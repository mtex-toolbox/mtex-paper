%% Weather Data: HAMLS versus LSQR
% This script computes the weather example of Section 3.1 of the paper. The
% data are the daily mean temperatures of 21 March 2018 at 6501 stations of
% the Global Historical Climatology Network. We compute the MLS reconstruction
% and the HAMLS approximation of bandwidth 256 and compare them with LSQR
% approximations that are given multiples of the HAMLS runtime.
%
% All output is written to the folder |results|:
%
%  WeatherData.png                          - the data (Figure 1a)
%  MLSApproximationWeatherData.png          - MLS reconstruction (Figure 1b)
%  HarmMLSApproximationWeatherData.png      - HAMLS approximation (Figure 1c)
%  LSQRApproximationWeatherData_<f>xT.png   - LSQR within f times the HAMLS
%                                             runtime (Figure 1d-g)
%  ColorbarWeather.png                      - color bar of Figure 1
%  IterationStudy_reg_*.txt                 - data of Figure 2
%  LSQRPanels.txt                           - runtimes and iteration counts
%
% The script has to be run from the folder |weather|. The data files are
% prepared by |prepareWeatherData.m|.
%
% The regularization parameter of LSQR cannot be optimized here, since no
% reference function is known. We use $\alpha = 10^{-7}$, chosen by
% inspection at bandwidth 128, and scale it to bandwidth 256 by $L^{-2s}$.

clear
close all

bw = 256;                    % bandwidth of the harmonic approximations
reg = 6.3e-9;                % regularization parameter of LSQR
sobolevIndex = 2;            % Sobolev index of the regularization
panelSize = [10 10 700 700]; % size of the sphere plots

dataDir = 'data';
outDir = 'results';
if ~isfolder(outDir), mkdir(outDir); end

%% Import the data
% Stations with identical positions are merged in the same way as
% |S2FunHarmonic.interpolate| does it, such that all quantities computed at
% the stations refer to the nodes the solvers actually use.

D = load(fullfile(dataDir,'Weather_Data.mat'));  % nodes, val
C = load(fullfile(dataDir,'coastLines.mat'));    % coastLines

assert(numel(D.nodes) == 6501,'expected 6501 stations');
assert(~any(isnan(D.val(:))),'the data contain NaN values');

[nodes,values] = uniqueStations(D.nodes(:),D.val(:));

% north pole up, prime meridian out of the screen
rF = specimenFrame.specimen;
rF.how2plot = 'z↑→y';
nodes.frame = rF;
coastLines = C.coastLines(:);
coastLines.frame = rF;

% all plots share this color range
colorRange = [min(values), max(values)];

fprintf('%i stations (%i after merging), %.1f to %.1f degrees Celsius\n', ...
  numel(D.nodes),numel(nodes),colorRange(1),colorRange(2));

% the degree at which the Sobolev penalty reaches the data term
lEff = reg^(-1/(2*sobolevIndex));
fprintf('LSQR: the penalty dominates from degree %.0f of %i on\n',lEff,bw);

%% Plot the data (Figure 1a)

plotStations(nodes,values,coastLines,colorRange);
savePng(fullfile(outDir,'WeatherData.png'),panelSize);

%% MLS reconstruction and HAMLS approximation (Figure 1b and 1c)
% The runtime |tHAMLS| covers Steps 1 to 3 of Algorithm 1 and is the time
% budget the LSQR approximations below refer to. The Voronoi areas cost only
% about 0.03 seconds here and are included in the runtimes of both methods.
%
% The MLS reconstruction uses tangent monomials of degree 4 and the target
% neighbor count $n = 4 \dim = 60$. The thresholds of the regularization
% (Section 2.3) are fixed instead of calibrated automatically. The warning of
% |S2FunMLS/eval| that |candidateFactor| should be raised can be ignored.

nRep = 10;      % the median over nRep runs is reported
mlsDegree = 4;  % polynomial degree of MLS
oF = 4;         % oversampling factor of MLS

mlsOpts = {'degree',mlsDegree,'oF',oF,'tangent',true,'weight','wendland', ...
  'use_vor_weights',true,'regularize',true, ...
  'mincond',10,'maxcond',30,'targetcond',10};

% warm up
S2FunMLS(nodes,values,mlsOpts{:});

[tMLS,mls] = medianTime(@() S2FunMLS(nodes,values,mlsOpts{:}),nRep);
[tHAMLS23,hamls] = medianTime(@() S2FunHarmonic(mls,'bandwidth',bw),nRep);
tHAMLS = tMLS + tHAMLS23;

assert(mls.oF == oF,'unexpected oversampling factor %g',mls.oF);

fprintf('MLS: %.2f s, fill distance %.1f degrees\n', ...
  tMLS,mls.fill_distance/degree);
fprintf('HAMLS: %.2f s for Steps 2 and 3, %.2f s in total\n',tHAMLS23,tHAMLS);

plotApprox(mls,coastLines,colorRange);
savePng(fullfile(outDir,'MLSApproximationWeatherData.png'),panelSize);

plotApprox(hamls,coastLines,colorRange);
savePng(fullfile(outDir,'HarmMLSApproximationWeatherData.png'),panelSize);

%% Iteration study (Figure 2)
% We run LSQR for an increasing number of iterations and record the two terms
% of the regularized functional,
%
% $$\mathrm{misfit} = \|W^{1/2}(F\hat f - y)\|_2, \qquad
%   \mathrm{penalty} = \sqrt{\alpha}\,\|\hat f\|_R.$$
%
% Since |lsqr| can neither be warm started nor be stopped by a callback, every
% iteration count requires a run of its own. This takes about 40 minutes.

regSweep = [1e-7 1e-8 1e-9 1e-10];
budgets = unique(round(logspace(1,log10(8000),12)));

nR = numel(regSweep);
nK = numel(budgets);

sweep = struct('reg',regSweep,'budget',budgets(:), ...
  'time',nan(nK,nR),'iter',nan(nK,nR),'flag',nan(nK,nR), ...
  'misfit',nan(nK,nR),'penalty',nan(nK,nR),'rms',nan(nK,nR));

% warm up
lsqrRun(nodes,values,reg,2,bw,sobolevIndex);

for i = 1:nR
  for k = 1:nK
    out = lsqrRun(nodes,values,regSweep(i),budgets(k),bw,sobolevIndex);

    sweep.time(k,i) = out.time;
    sweep.iter(k,i) = out.iter;
    sweep.flag(k,i) = out.flag;
    sweep.misfit(k,i) = out.misfit;
    sweep.penalty(k,i) = out.penalty;
    sweep.rms(k,i) = out.rms;

    fprintf(['alpha = %-8g %5i iterations, %6.1f s, misfit %.4g, ', ...
      'penalty %.4g, rms %.2f\n'],regSweep(i),out.iter,out.time, ...
      out.misfit,out.penalty,out.rms);
  end
end

save(fullfile(outDir,'IterationStudy.mat'),'sweep','bw','sobolevIndex');

% One table per regularization parameter. The file names write the minus as
% 'm' and the decimal point as 'p', e.g. IterationStudy_reg_6p3em09.txt.
% solutionNorm is the penalty without the factor sqrt(alpha).
for i = 1:nR
  tag = strrep(strrep(sprintf('%g',regSweep(i)),'-','m'),'.','p');
  writeTable(fullfile(outDir,sprintf('IterationStudy_reg_%s.txt',tag)), ...
    {'iterations','misfit','penalty','solutionNorm','rms','time'}, ...
    [sweep.iter(:,i), sweep.misfit(:,i), sweep.penalty(:,i), ...
    sweep.penalty(:,i)/sqrt(regSweep(i)), sweep.rms(:,i), sweep.time(:,i)]);
end

% a quick look at the result; Figure 2 is drawn from the tables
figure
loglog(sweep.iter,sweep.misfit,'-o')
hold on
set(gca,'ColorOrderIndex',1)
loglog(sweep.iter,sweep.penalty,'--')
hold off
grid on
xlabel('LSQR iterations')
ylabel('misfit (solid) and penalty (dashed)')
legend(arrayfun(@(a) sprintf('\\alpha = %g',a),regSweep, ...
  'UniformOutput',false),'Location','best')

%% LSQR within multiples of the HAMLS runtime (Figure 1d-g)
% LSQR is given a time budget of |timeFactors(j)| times the HAMLS runtime.
% The number of iterations that fits into this budget is predicted by the
% affine model $t(K) = t_0 + cK$. As this model is only approximate, a run
% that misses its budget by more than |timeTol| is repeated with a corrected
% number of iterations. The paper shows the factors 1, 5, 10 and 30.

timeFactors = [1 5 10 15 20 25 30 35 40 50 100];
calibIt = 200;   % iterations of the calibration run
timeTol = 0.05;  % accepted relative deviation from the budget
maxRefine = 3;   % maximum number of corrections

% The runtime of HAMLS is measured again, since the iteration study took
% almost an hour and the machine may have warmed up in the meantime.
tHAMLS23new = medianTime(@() S2FunHarmonic(mls,'bandwidth',bw),nRep);
if abs(tHAMLS23new - tHAMLS23) > 0.1 * tHAMLS23
  warning('The HAMLS runtime changed from %.2f s to %.2f s.', ...
    tHAMLS23,tHAMLS23new);
end
tHAMLS23 = tHAMLS23new;
tHAMLS = tMLS + tHAMLS23;

% calibrate the time model at 0 and calibIt iterations
calib0 = lsqrRun(nodes,values,reg,0,bw,sobolevIndex);
calibK = lsqrRun(nodes,values,reg,calibIt,bw,sobolevIndex);
t0 = calib0.time;
cPerIt = (calibK.time - t0) / calibIt;
assert(cPerIt > 0,'the LSQR runtime does not grow with the iterations');
fprintf('LSQR: %.2f s setup and %.1f ms per iteration\n',t0,1e3*cPerIt);

panel = struct('name',{},'iter',{},'time',{},'misfit',{},'rms',{}, ...
  'range',{},'f',{});

for j = 1:numel(timeFactors)
  target = timeFactors(j) * tHAMLS;
  K = budgetFor(target,t0,cPerIt);

  for correction = 0:maxRefine
    out = lsqrRun(nodes,values,reg,K,bw,sobolevIndex);
    fprintf('%gx: %.1f s budget, %5i iterations, %.1f s (%+.0f %%)\n', ...
      timeFactors(j),target,out.iter,out.time,100*(out.time-target)/target);

    if abs(out.time - target) <= timeTol * target || out.time <= t0, break, end

    % correct the model by the cost per iteration of this run
    K = budgetFor(target,t0,(out.time - t0) / K);
  end

  panel(j) = summarizePanel(sprintf('LSQR_%gxT',timeFactors(j)),out);

  plotApprox(out.f,coastLines,colorRange);
  savePng(fullfile(outDir, ...
    sprintf('LSQRApproximationWeatherData_%gxT.png',timeFactors(j))),panelSize);
end

%% LSQR at the iteration counts of the paper (optional)
% The iteration counts of the previous section depend on the machine. With
% |paperCounts = true| the LSQR panels are redrawn at the iteration counts
% stored in |results/LSQRPanels.txt|, i.e., the counts of Figure 1d-g. This
% has to be done before the next section overwrites that file.

paperCounts = false;

if paperCounts
  paperPanels = readtable(fullfile(outDir,'LSQRPanels.txt')); %#ok<UNRCH>
  isPanel = ~cellfun(@isempty,regexp(paperPanels.panel,'^LSQR_[\d.]+xT$','once'));

  for row = find(isPanel)'
    name = paperPanels.panel{row};
    factor = sscanf(name,'LSQR_%gxT');

    out = lsqrRun(nodes,values,reg,paperPanels.iterations(row),bw,sobolevIndex);
    fprintf('%s: %i iterations, %.1f s\n',name,out.iter,out.time);

    plotApprox(out.f,coastLines,colorRange);
    savePng(fullfile(outDir, ...
      sprintf('LSQRApproximationWeatherData_%gxT.png',factor)),panelSize);
  end
end

%% Color bar
% The color bar is cut out of this plot in the paper, so the size of the
% figure should not be changed.

figure
plot(hamls,'upper','nolabel');
colormap(WhiteJetColorMap);
setColorRange(colorRange);
mtexColorbar;
savePng(fullfile(outDir,'ColorbarWeather.png'),[10 10 1000 1200]);

%% Summary
% We print the numbers quoted in Section 3.1 and store them in
% |LSQRPanels.txt|. The range of every approximation is determined on a grid
% of resolution one degree.

w = voronoiWeights(nodes);
gEval = equispacedS2Grid('resolution',1*degree);
rangeOf = @(f) [min(f.eval(gEval)), max(f.eval(gEval))];
rmsAt = @(f) sqrt(mean((f.eval(nodes) - values).^2));

rowMLS = summarizePanel('MLS',struct('iter',NaN,'time',tMLS, ...
  'misfit',norm(sqrt(w) .* (mls.eval(nodes) - values)), ...
  'rms',rmsAt(mls),'f',mls));
rowHAMLS = summarizePanel('HAMLS',struct('iter',NaN,'time',tHAMLS, ...
  'misfit',norm(sqrt(w) .* (hamls.eval(nodes) - values)), ...
  'rms',rmsAt(hamls),'f',hamls));

rows = [rowMLS, rowHAMLS, panel];
for k = 1:numel(rows)
  rows(k).range = rangeOf(rows(k).f);
end

fprintf('\nSection 3.1\n');
fprintf('stations         %i\n',numel(nodes));
fprintf('fill distance    %.1f degrees\n',mls.fill_distance/degree);
fprintf('bandwidth        %i\n',bw);
fprintf('data range       [%.1f, %.1f]\n\n',colorRange);
fprintf('%-12s %8s %8s %10s %8s %18s\n', ...
  '','iter','seconds','misfit','rms','range');
for k = 1:numel(rows)
  fprintf('%-12s %8g %8.1f %10.4g %8.2f  [%6.1f, %6.1f]\n', ...
    rows(k).name,rows(k).iter,rows(k).time,rows(k).misfit,rows(k).rms, ...
    rows(k).range);
end

writeTable(fullfile(outDir,'LSQRPanels.txt'), ...
  {'panel','iterations','seconds','misfit','rms','minimum','maximum'}, ...
  [reshape([rows.iter],[],1), reshape([rows.time],[],1), ...
  reshape([rows.misfit],[],1), reshape([rows.rms],[],1), ...
  reshape(cellfun(@(r) r(1),{rows.range}),[],1), ...
  reshape(cellfun(@(r) r(2),{rows.range}),[],1)],{rows.name});

paperNumbers = struct('stations',numel(nodes), ...
  'fillDeg',mls.fill_distance/degree,'bandwidth',bw, ...
  'mlsTime',tMLS,'hamlsTime23',tHAMLS23,'hamlsTime',tHAMLS, ...
  'lsqrPerIt',cPerIt,'reg',reg,'lEff',lEff);

% assigned separately, since struct() would turn a struct array value into a
% struct array
paperNumbers.panels = rows;

%% Helper functions

function out = lsqrRun(nodes,values,alpha,K,L,s)
% LSQR with K iterations and the two terms of the regularized functional
%
% Syntax
%   out = lsqrRun(nodes,values,alpha,K,L,s)
%
% Input
%  nodes, values - scattered data
%  alpha - regularization parameter
%  K     - number of iterations
%  L     - bandwidth
%  s     - Sobolev index
%
% Output
%  out - struct with fields f, time, flag, iter, misfit, penalty, rms
%
% The tolerance is chosen out of reach, such that LSQR performs exactly K
% iterations. Otherwise a note is printed.

ws = warning('off','lsqr:itermax');
tic
[f,p] = S2FunHarmonic.interpolate(nodes,values,'bandwidth',L, ...
  'regularization',alpha,'SobolevIndex',s,'maxit',K,'tol',1e-15);
out.time = toc;
warning(ws);

out.f = f;
out.flag = p{1};
out.iter = p{3};

if out.flag ~= 1 && K > 0
  fprintf('  note: lsqr stopped with flag %i after %i of %i iterations\n', ...
    out.flag,out.iter,K);
end

w = voronoiWeights(nodes);
What = sobolevWeights(L,s);

res = f.eval(nodes) - values;
out.misfit = norm(sqrt(w) .* res);
out.penalty = sqrt(alpha) * norm(sqrt(What(1:numel(f.fhat))) .* f.fhat(:));
out.rms = sqrt(mean(res.^2));

% The two terms are the two blocks of the residual that lsqr minimizes, so
% their Pythagorean sum has to equal the residual norm reported by lsqr. This
% checks the weights and alpha. Since lsqr updates this norm recursively, the
% tolerance is loose.
resvec = p{4}{1};
if K > 0 && abs(hypot(out.misfit,out.penalty) - resvec(end)) > 1e-3 * resvec(end)
  warning('misfit and penalty do not match the residual of lsqr (%.6g vs %.6g)', ...
    hypot(out.misfit,out.penalty),resvec(end));
end

end


function K = budgetFor(target,t0,cPerIt)
% number of iterations that fit into target seconds, at least one

K = max(1,round((target - t0) / cPerIt));

end


function w = voronoiWeights(nodes)
% normalized Voronoi areas, the weights used by S2FunHarmonic.interpolate
%
% As the weights sum up to one, the misfit is a weighted root mean square
% error in degrees Celsius.

w = calcVoronoiArea(nodes) / 4 / pi;
w = w(:);
assert(abs(sum(w) - 1) < 1e-9,'the Voronoi weights do not sum up to one');

end


function What = sobolevWeights(L,s)
% Sobolev weights (1 + l(l+1))^s as used by S2FunHarmonic.interpolate

What = repelem((1 + (0:L) .* ((0:L) + 1)).^s, 1:2:(2*L+1))';

end


function [nodes,values] = uniqueStations(nodes,values)
% merge stations at identical positions and average their values

[nodes,~,ind] = unique(nodes(:));
values = accumarray(ind,values(:),[]) ./ accumarray(ind,1);

end


function [t,out] = medianTime(fun,nRep)
% median runtime of nRep calls of fun and the result of the last call

ts = nan(1,nRep);
for r = 1:nRep
  tic
  out = fun();
  ts(r) = toc;
end
t = median(ts,'omitnan');

end


function row = summarizePanel(name,out)
% one row of the table LSQRPanels.txt

row = struct('name',name,'iter',out.iter,'time',out.time, ...
  'misfit',out.misfit,'rms',out.rms,'range',[NaN NaN],'f',out.f);

end


function plotStations(nodes,values,coastLines,colorRange)
% scatter plot of the data together with the coast lines

figure
plot(nodes,values,'upper','nolabel','markersize',5);
setColorRange(colorRange);
colormap(WhiteJetColorMap);
hold on
plot(coastLines,'upper','nolabel','markersize',1, ...
  'markeredgecolor','k','markerfacecolor','k');
hold off

end


function plotApprox(f,coastLines,colorRange)
% plot of an approximation together with the coast lines

figure
plot(f,'upper','nolabel','colorRange',colorRange);
colormap(WhiteJetColorMap);
hold on
plot(coastLines,'upper','nolabel','markersize',1, ...
  'markeredgecolor','k','markerfacecolor','k');
hold off

end


function savePng(fname,figSize)
% resize the current figure and save it as png
%
% MTEX keeps its axes at a fixed size in pixels and adjusts them in the
% ResizeFcn of the figure. Without a display MATLAB does not call this
% function, so we call it explicitly.

set(gcf,'position',figSize);
resizeFcn = get(gcf,'ResizeFcn');
if ~isempty(resizeFcn), feval(resizeFcn,gcf,[]); end

saveas(gcf,fname,'png');

end


function writeTable(fname,header,data,rowNames)
% write a tab separated table with a single header line
%
% Syntax
%   writeTable(fname,header,data)
%   writeTable(fname,header,data,rowNames)

assert(numel(header) == size(data,2) + (nargin > 3), ...
  '%s: %i column names for %i columns',fname,numel(header),size(data,2));

fid = fopen(fname,'w');
fprintf(fid,'%s',header{1});
fprintf(fid,'\t%s',header{2:end});
fprintf(fid,'\n');
if nargin > 3
  for k = 1:size(data,1)
    fprintf(fid,'%s',rowNames{k});
    fprintf(fid,'\t%.10g',data(k,:));
    fprintf(fid,'\n');
  end
end
fclose(fid);

if nargin <= 3
  writematrix(data,fname,'Delimiter','tab','WriteMode','append');
end

end
