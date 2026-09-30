%% Synthetic Example: Error over the Bandwidth
% This script computes the synthetic example of Section 3.2 of the paper. A
% test function of bandwidth 512 is sampled at $N = 10^4$ strongly
% non-uniform nodes. From these data we compute HAMLS approximations of
% polynomial degree 1 to 4 and LSQR approximations with an error minimizing
% regularization parameter, for all bandwidths from 4 to 128.
%
% All output is written to the folder |results/errorVsBandwidth|:
%
%  ToyFun.png, ToyData.png, ColorbarToy.eps - Figure 3
%  ApproximationError.txt                   - Figure 4
%  ApproximationTime.txt                    - runtimes quoted in Section 3.2
%  ToyExampleData.mat                       - the data set
%  RegParameters.mat                        - the regularization parameters
%  SetupTime.mat, LSQRError.mat, MLSError.mat - individual measurements
%
% The data set and the regularization parameters are used by
% |equalTimeError.m| and |stepRuntimes.m|, so this script has to be run
% first. It has to be run from the folder |synthetic| and takes several
% hours, most of it for the search for the regularization parameters.
%
% The runtimes include the computation of the Voronoi areas of the nodes,
% which both methods need. The runtime of HAMLS is the construction of the
% @S2FunMLS plus |S2FunHarmonic(mls,'bandwidth',L)|, i.e., Steps 1 to 3 of
% Algorithm 1.

clear
close all

dataDir = 'data';
outDir = fullfile('results','errorVsBandwidth');
if ~isfolder(outDir), mkdir(outDir); end

bw = 4:128;       % bandwidths
degrees = 1:4;    % polynomial degrees of MLS
nRep = 5;         % the median over nRep runs is reported

nB = numel(bw);
nD = numel(degrees);

% options of every LSQR run
lsqrOpts = {'SobolevIndex',2,'maxit',10000,'tol',1e-6};

%% The test function and the nodes
% The nodes are drawn from a smoothed version of the test function, which
% |discreteSample| normalizes to a probability density. A fraction of 4
% percent of the nodes is drawn uniformly, such that the gaps do not become
% arbitrarily large. |stepRuntimes.m| draws the same nodes again and checks
% them against |ToyExampleData.mat|.

nNodes = 1e4;
uniformFraction = 0.04;

F = load(fullfile(dataDir,'ToyExampleFun.mat'));  % fun
fun = F.fun * 8 * pi^2 / sum(F.fun);

% the density is the smoothed function, the data are the shifted function
g = smooth(fun,'halfwidth',5*degree);
fun = fun + 300;

nUniform = round(uniformFraction * nNodes);
rng(3)
nodes = [discreteSample(g,nNodes - nUniform); vector3d.rand(nUniform)];
values = fun.eval(nodes);

save(fullfile(outDir,'ToyExampleData'),'fun','nodes','values');

%% MLS approximations
% The construction of the MLS approximation does not depend on the bandwidth.
% Its runtime, which includes the Voronoi areas, is added to every HAMLS
% runtime below. The Voronoi areas are timed separately as well.

mlsOpts = {'oF',4,'tangent',true,'regularize',true,'weight','auto', ...
  'use_smooth_delta',true,'use_vor_weights',true};

% warm up, as the first call carries JIT and allocation overhead
calcVoronoiArea(nodes);
[tVormed,~,tVor] = medianTime(@() calcVoronoiArea(nodes),nRep);

S2FunMLS(nodes,values,'degree',degrees(1),mlsOpts{:});

mls = cell(nD,1);
tSetup = nan(nD,nRep);
tSetupmed = nan(nD,1);
for j = 1:nD
  [tSetupmed(j),mls{j},tSetup(j,:)] = medianTime( ...
    @() S2FunMLS(nodes,values,'degree',degrees(j),mlsOpts{:}),nRep);
end

fprintf('Voronoi areas: %.4f s\n',tVormed);
fprintf('MLS construction: %s s\n',mat2str(round(tSetupmed',4)));
fprintf('fill distance: %.2f degrees\n',mls{1}.fill_distance/degree);

save(fullfile(outDir,'SetupTime'),'degrees','nRep','tVor','tSetup', ...
  'tVormed','tSetupmed');

%% L2-projection
% The $L_2$-projection of the test function onto the functions of bandwidth
% $L$ is the truncation of its Fourier coefficients. It is the best
% approximation of bandwidth $L$ and serves as a lower bound. Note that the
% bandwidth has to be set as a property, since the constructor ignores the
% option |'bandwidth'| for an @S2FunHarmonic.

assert(fun.bandwidth >= bw(end),'the test function has bandwidth %i',fun.bandwidth);

ErrorProj = zeros(nB,1);
for n = 1:nB
  fL = fun;
  fL.bandwidth = bw(n);
  ErrorProj(n) = relError(fL,fun);
end

%% LSQR: the optimal regularization parameter
% For every bandwidth we choose the regularization parameter that minimizes
% the $L_2$-error. This search is by far the most expensive part of the
% script. Its runtime is not counted as runtime of LSQR.

nRefine = 4;  % refinement steps of the search

Reg = nan(nB,1);
ws = warning('off','all');
for n = 1:nB
  [Reg(n),err,nSolves] = findBestRegularization(nodes,values,fun, ...
    bw(n),lsqrOpts,nRefine);
  fprintf('bandwidth %3i: reg = %.3g, error = %.4g (%i solves)\n', ...
    bw(n),Reg(n),err,nSolves);
end
warning(ws);

save(fullfile(outDir,'RegParameters'),'bw','Reg');

%% LSQR: error and runtime

tLSQR = nan(nB,nRep);
tLSQRmed = nan(nB,1);
ErrorLSQR = nan(nB,1);
lsqrParameters = cell(nB,1);

for n = 1:nB
  fprintf('LSQR, bandwidth %i\n',bw(n));
  [tLSQRmed(n),fit,tLSQR(n,:)] = medianTime( ...
    @() lsqrFit(nodes,values,bw(n),Reg(n),lsqrOpts),nRep);

  ErrorLSQR(n) = relError(fit.f,fun);
  lsqrParameters{n} = fit.par;
end

save(fullfile(outDir,'LSQRError'),'bw','Reg','ErrorLSQR','tLSQR', ...
  'tLSQRmed','lsqrParameters');

%% HAMLS: error and runtime
% Here we time only Steps 2 and 3. The construction of the MLS approximation
% is added when the tables are written.

tConv = nan(nB,nD,nRep);
tConvmed = nan(nB,nD);
ErrorMLS = nan(nB,nD);

for n = 1:nB
  fprintf('HAMLS, bandwidth %i\n',bw(n));
  for j = 1:nD
    [tConvmed(n,j),harmls,tConv(n,j,:)] = medianTime( ...
      @() S2FunHarmonic(mls{j},'bandwidth',bw(n)),nRep);

    ErrorMLS(n,j) = relError(harmls,fun);
  end
end

save(fullfile(outDir,'MLSError'),'bw','degrees','ErrorMLS','tConv','tConvmed');

%% Plot the test function and the nodes (Figure 3)

colorRange = [300 410];

figure(1)
plot(fun,'upper','nolabel','colorRange',colorRange);
colormap(WhiteJetColorMap);
fcw;
set(gcf,'position',[10 10 800 800]);
exportgraphics(gcf,fullfile(outDir,'ToyFun.png'),'Resolution',300);

figure(2)
plot(nodes,values,'upper','nolabel','markersize',5,'all','MarkerEdgeColor','k');
colormap(WhiteJetColorMap);
fcw;
set(gcf,'position',[10 10 500 500]);
exportgraphics(gcf,fullfile(outDir,'ToyData.png'),'Resolution',300);

% The color bar is cut out of this plot in the paper, so the size of the
% figure should not be changed.
figure(3)
plot(fun,'upper','nolabel');
setColorRange(colorRange);
mtexColorbar;
set(gcf,'position',[10 10 600 500]);
saveas(gcf,fullfile(outDir,'ColorbarToy.eps'),'epsc');

%% Write the tables
% Figure 4 is drawn from |ApproximationError.txt| with pgfplots.

timeLSQR = tLSQRmed;
timeHAMLS = tConvmed + tSetupmed';

degNames = arrayfun(@(d) sprintf('MLS_Harm_deg%i',d),degrees,'UniformOutput',false);

writeTable(fullfile(outDir,'ApproximationError.txt'), ...
  [{'bandwidth','Projection','LSQR'},degNames],[bw', ErrorProj, ErrorLSQR, ErrorMLS]);

writeTable(fullfile(outDir,'ApproximationTime.txt'), ...
  [{'bandwidth','LSQR'},degNames],[bw', timeLSQR, timeHAMLS]);

%% Summary
% We print the numbers quoted in Section 3.2. The HAMLS errors stagnate from
% about the bandwidth on at which $(L+1)^2$ exceeds the number of nodes.
% Section 3.2 also quotes the fraction of the sphere on which the local MLS
% problems are not regularized. We measure it as area, by the weights of a
% quadrature grid, since the $L_2$-error averages over the area.

nLast = nB;                                 % L = 128
jPlot = find(ismember(degrees,[2 3 4]));    % the degrees shown in Figure 4

lsqrIterLast = lsqrParameters{nLast}{3};

iStag = find((bw + 1).^2 > nNodes,1);
if isempty(iStag), iStag = nB; end

S2G = quadratureS2Grid(bw(nLast),'GaussLegendre');
quadW = S2G.weights(:) / sum(S2G.weights(:));
regFree = nan(1,nD);
for j = jPlot
  [~,~,regInfo] = mls{j}.eval(S2G);
  regFree(j) = sum(quadW(~logical(regInfo.regularizationActive(:))));
end

fprintf('\nSection 3.2\n');
fprintf('nodes                   %i\n',numel(nodes));
fprintf('fill distance           %.1f degrees\n',mls{1}.fill_distance/degree);
fprintf('Voronoi areas           %.2f s\n',tVormed);
fprintf('MLS construction        %.2f to %.2f s\n', ...
  min(tSetupmed(jPlot)),max(tSetupmed(jPlot)));
fprintf('L = %i, LSQR           %.1f s, %i iterations\n', ...
  bw(nLast),timeLSQR(nLast),lsqrIterLast);
fprintf('L = %i, HAMLS          %.2f to %.2f s\n', ...
  bw(nLast),min(timeHAMLS(nLast,jPlot)),max(timeHAMLS(nLast,jPlot)));
fprintf('L = %i, errors         LSQR %.2e, HAMLS %.2e to %.2e\n', ...
  bw(nLast),ErrorLSQR(nLast),min(ErrorMLS(nLast,jPlot)),max(ErrorMLS(nLast,jPlot)));
fprintf('(L+1)^2 > N from        L = %i, HAMLS errors %.2e to %.2e\n', ...
  bw(iStag),min(ErrorMLS(iStag,jPlot)),max(ErrorMLS(iStag,jPlot)));
fprintf('not regularized         %.0f%% to %.0f%% of the sphere\n', ...
  100*min(regFree(jPlot)),100*max(regFree(jPlot)));

paperNumbers = struct('nodes',numel(nodes), ...
  'fillDeg',mls{1}.fill_distance/degree,'tVor',tVormed, ...
  'tSetupRange',[min(tSetupmed(jPlot)), max(tSetupmed(jPlot))], ...
  'lsqrTime128',timeLSQR(nLast),'lsqrIter128',lsqrIterLast, ...
  'hamlsRange128',[min(timeHAMLS(nLast,jPlot)), max(timeHAMLS(nLast,jPlot))], ...
  'stagnationL',bw(iStag), ...
  'stagnationErr',[min(ErrorMLS(iStag,jPlot)), max(ErrorMLS(iStag,jPlot))], ...
  'regFreeRange',[min(regFree(jPlot)), max(regFree(jPlot))]);

%% A quick look at the results

figure(11)
semilogy(bw,ErrorProj,'k-','LineWidth',1.5); hold on
semilogy(bw,ErrorLSQR,'r-','LineWidth',1.5);
semilogy(bw,ErrorMLS,'--');
hold off
xlabel('bandwidth L'); ylabel('relative L2-error'); grid on
legend([{'L2-projection','LSQR'}, ...
  arrayfun(@(d) sprintf('HAMLS (degree %i)',d),degrees,'UniformOutput',false)], ...
  'Location','southwest');

figure(12)
semilogy(bw,timeLSQR,'r-','LineWidth',1.5); hold on
semilogy(bw,timeHAMLS,'--');
hold off
xlabel('bandwidth L'); ylabel('runtime in seconds'); grid on
legend([{'LSQR'}, ...
  arrayfun(@(d) sprintf('HAMLS (degree %i)',d),degrees,'UniformOutput',false)], ...
  'Location','northwest');

%% Helper functions

function out = lsqrFit(nodes,values,L,alpha,opts)
% regularized least squares approximation by LSQR
%
% Output
%  out.f   - @S2FunHarmonic
%  out.par - lsqr diagnostics as returned by S2FunHarmonic.interpolate

[out.f,out.par] = S2FunHarmonic.interpolate(nodes,values, ...
  'bandwidth',L,'regularization',alpha,opts{:});

end


function [t,out,ts] = medianTime(fun,nRep)
% median runtime of nRep calls of fun, the result of the last call and all
% runtimes

ts = nan(1,nRep);
for r = 1:nRep
  tic
  out = fun();
  ts(r) = toc;
end
t = median(ts,'omitnan');

end


function e = relError(f,ref)
% relative L2-error of two S2FunHarmonic, computed from their coefficients
%
% We do not use norm(f - ref) here. The difference of two S2FunHarmonic is
% truncated to the coefficients larger than 1e-8 times the largest one, which
% for the test function of this script means bandwidth 112 instead of 512.

a = f.fhat(:);
b = ref.fhat(:);
n = max(numel(a),numel(b));
a(end+1:n,1) = 0;
b(end+1:n,1) = 0;
e = norm(a - b) / norm(b);

end


function [reg,err,nSolves] = findBestRegularization(nodes,values,fun,L,opts,nRefine)
% regularization parameter that minimizes the L2-error at bandwidth L
%
% Syntax
%   [reg,err,nSolves] = findBestRegularization(nodes,values,fun,L,opts,nRefine)
%
% Output
%  reg     - best regularization parameter
%  err     - relative L2-error at reg
%  nSolves - number of LSQR runs
%
% A coarse logarithmic grid locates the minimum, then the bracket around it
% is bisected nRefine times on the logarithmic scale.

nSolves = 0;

% coarse grid, decreasing; since the error is U-shaped in the parameter we
% stop once we are well past the minimum
r = 10.^-(0:16);
E = nan(size(r));
for id = 1:numel(r)
  E(id) = errorFor(r(id));
  if id > 4 && E(id) > 1.3 * min(E), break; end
end
nEval = find(~isnan(E),1,'last');

[~,id] = min(E);

% warnings are switched off by the caller, so this is printed
if id == 1 || id == nEval
  fprintf('  note: the optimum is at the boundary of the coarse grid (bw = %i, reg = %.3g)\n', ...
    L,r(id));
end

% bracket the minimum on the log10 scale, clamped to the coarse grid
id = min(max(id,2),nEval-1);
hi = log10(r(id-1));
m = log10(r(id));
lo = log10(r(id+1));
Em = E(id);

for k = 1:nRefine
  c1 = (lo + m) / 2;
  c2 = (m + hi) / 2;
  E1 = errorFor(10^c1);
  E2 = errorFor(10^c2);
  [Em,w] = min([E1, Em, E2]);
  switch w
    case 1, hi = m;  m = c1;
    case 2, lo = c1; hi = c2;
    case 3, lo = m;  m = c2;
  end
end

reg = 10^m;
err = Em;

  function e = errorFor(alpha)
    e = relError(lsqrFit(nodes,values,L,alpha,opts).f,fun);
    nSolves = nSolves + 1;
  end

end


function writeTable(fname,header,data)
% write a tab separated table with a single header line

assert(numel(header) == size(data,2), ...
  '%s: %i column names for %i columns',fname,numel(header),size(data,2));

fid = fopen(fname,'w');
fprintf(fid,'%s',header{1});
fprintf(fid,'\t%s',header{2:end});
fprintf(fid,'\n');
fclose(fid);

writematrix(data,fname,'Delimiter','tab','WriteMode','append');

end
