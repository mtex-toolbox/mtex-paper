%% Synthetic Example: LSQR within Multiples of the HAMLS Runtime
% This script computes Figure 5 of the paper. For every bandwidth from 4 to
% 128 we compare the error of the HAMLS approximation of polynomial degree 4
% with the error LSQR achieves when it is given |timeFactors| times the
% runtime of HAMLS.
%
% All output is written to the folder |results/equalTime|:
%
%  LSQRErrorTimeFactors.txt - errors, drawn in Figure 5
%  LSQRIterationBudgets.txt - HAMLS runtimes and LSQR iteration budgets
%  EqualTimeHAMLS.mat, EqualTimeLSQRCalibration.mat,
%  EqualTimeErrorTimeFactors.mat - individual measurements
%
% The script reads the data set and the regularization parameters written by
% |errorVsBandwidth.m|. It has to be run from the folder |synthetic| and takes
% about two days; an estimate is printed before the expensive part starts.
%
% Both methods compute the Voronoi areas of the nodes. Since this step does
% not depend on the bandwidth, it is timed once and subtracted from both
% runtimes. Hence we compare
%
% $$t_{\mathrm{HAMLS}}(L) = t_{\mathrm{MLS}} + t_{\mathrm{Steps\,2,3}}(L)
%   - t_{\mathrm{Voronoi}}$$
%
% with the affine model $t_{\mathrm{LSQR}}(K) = t_0(L) - t_{\mathrm{Voronoi}}
% + c(L) K$ for $K$ iterations, which is calibrated for every bandwidth at
% $K = 0$ and $K = 100$. The iteration budget of the factor $f$ is the largest
% $K$ with $t_{\mathrm{LSQR}}(K) \le f\, t_{\mathrm{HAMLS}}(L)$. If not a single
% iteration fits into the budget, LSQR returns its initial guess and the
% relative error is 1.
%
% Stopping early acts as a regularization itself, so the optimal
% regularization parameter depends on the budget. We therefore search it for
% every bandwidth and factor, starting from the parameter found for the
% previous factor, or the previous bandwidth, or the fully iterated problem.

clear
close all

toyDir = fullfile('results','errorVsBandwidth');
outDir = fullfile('results','equalTime');
if ~isfolder(outDir), mkdir(outDir); end

bw = 4:128;                                    % bandwidths
mlsDegree = 4;                                 % polynomial degree of MLS
timeFactors = [1 5 10 15 20 25 30 35 40 50 100];  % time budgets of LSQR
nRep = 3;                                      % the median over nRep runs is used

nB = numel(bw);
nF = numel(timeFactors);

% options of every LSQR run; the tolerance is out of reach, such that only
% the iteration budget stops LSQR
lsqrOpts = {'SobolevIndex',2,'tol',1e-15};

assert(issorted(timeFactors) && all(timeFactors > 0), ...
  'timeFactors has to be positive and increasing');

% The timings are taken on a fixed set of CPUs (Linux only). '0-15' with 16
% threads uses all hardware threads of an 8-core CPU with simultaneous
% multithreading, like the one of Section 3. On a CPU with performance and
% efficiency cores choose the performance cores, since the cost of an LSQR
% iteration determines every budget. Set pinCores = '' to leave the choice
% to the operating system.
pinCores = '0-15';
nThreads = 16;
if ~isempty(pinCores) && isunix
  system(sprintf('taskset -a -pc %s %d > /dev/null 2>&1',pinCores,feature('getpid')));
  maxNumCompThreads(nThreads);
end

%% Import the data set of errorVsBandwidth.m

D = load(fullfile(toyDir,'ToyExampleData.mat'));  % fun, nodes, values
fun = D.fun;
nodes = D.nodes;
values = D.values;

assert(fun.bandwidth >= bw(end),'the test function has bandwidth %i',fun.bandwidth);

%% MLS approximation and Voronoi areas
% The options are those of |errorVsBandwidth.m|, such that the runtimes
% measured here belong to the approximations of Figure 4.

mlsOpts = {'oF',4,'tangent',true,'regularize',true,'weight','auto', ...
  'use_smooth_delta',true,'use_vor_weights',true};

% warm up, as the first call carries JIT and allocation overhead
calcVoronoiArea(nodes);
[tVormed,~,tVor] = medianTime(@() calcVoronoiArea(nodes),nRep);

S2FunMLS(nodes,values,'degree',mlsDegree,mlsOpts{:});
[tSetupmed,mls,tSetup] = medianTime( ...
  @() S2FunMLS(nodes,values,'degree',mlsDegree,mlsOpts{:}),nRep);

fprintf('Voronoi areas: %.4f s\n',tVormed);
fprintf('MLS construction: %.4f s\n',tSetupmed);

%% HAMLS: error and runtime
% |S2FunHarmonic(mls,'bandwidth',L)| performs Steps 2 and 3, the evaluation
% of the MLS approximation on the quadrature grid and the adjoint spherical
% Fourier transform.

tConv = nan(nB,nRep);
tConvmed = nan(nB,1);
ErrorMLS = nan(nB,1);

for n = 1:nB
  fprintf('HAMLS, bandwidth %i\n',bw(n));
  [tConvmed(n),harmls,tConv(n,:)] = medianTime( ...
    @() S2FunHarmonic(mls,'bandwidth',bw(n)),nRep);

  ErrorMLS(n) = relError(harmls,fun);
end

save(fullfile(outDir,'EqualTimeHAMLS'),'bw','mlsDegree','nRep','ErrorMLS', ...
  'tConv','tConvmed','tSetup','tSetupmed','tVor','tVormed');

% runtime of HAMLS without the Voronoi areas
tHAMLS = tConvmed + tSetupmed - tVormed;
assert(all(tHAMLS > 0),'nonpositive HAMLS runtime');

%% L2-projection
% The truncation of the Fourier coefficients of the test function is a lower
% bound for every approximation of bandwidth $L$.

ErrorProj = zeros(nB,1);
for n = 1:nB
  fL = fun;
  fL.bandwidth = bw(n);
  ErrorProj(n) = relError(fL,fun);
end

%% Time model of LSQR and iteration budgets
% One LSQR iteration costs a nonequispaced spherical Fourier transform and
% its adjoint. Its cost depends strongly, and not monotonically, on the
% bandwidth, so the time model is calibrated for every bandwidth.

calibIt = 100;    % iterations of the calibration run
maxItCap = 1e5;   % upper bound for the budgets
regCalib = 1e-7;  % the regularization does not change the cost of an iteration

tLSQR0 = nan(nB,nRep);    % LSQR with zero iterations
tLSQRcal = nan(nB,nRep);  % LSQR with calibIt iterations

ws = warning('off','all');
for n = 1:nB
  fprintf('LSQR calibration, bandwidth %i\n',bw(n));
  for r = 1:nRep
    [~,tLSQR0(n,r)] = lsqrFit(nodes,values,bw(n),regCalib,0,lsqrOpts);
    [~,tLSQRcal(n,r)] = lsqrFit(nodes,values,bw(n),regCalib,calibIt,lsqrOpts);
  end
end
warning(ws);

save(fullfile(outDir,'EqualTimeLSQRCalibration'),'bw','nRep','calibIt', ...
  'tLSQR0','tLSQRcal');

t0 = median(tLSQR0,2,'omitnan');
tCal = median(tLSQRcal,2,'omitnan');
slope = (tCal - t0) / calibIt;  % seconds per iteration
t0net = t0 - tVormed;           % setup without the Voronoi areas

% For small bandwidths the difference of t0 and tCal may be below the timing
% noise, and the slope may even be negative. We drop these bandwidths.
bad = slope <= 0;
if any(bad)
  warning('The time model failed for the bandwidths %s; increase nRep or calibIt.', ...
    mat2str(bw(bad)));
end

iterBudget = nan(nB,nF);
for k = 1:nF
  iterBudget(:,k) = floor(max(0,(timeFactors(k) * tHAMLS - t0net) ./ slope));
end
iterBudget(bad,:) = NaN;

if any(iterBudget(:) > maxItCap)
  warning('%i budgets are cut to %g iterations.',sum(iterBudget(:) > maxItCap),maxItCap);
  iterBudget = min(iterBudget,maxItCap);
end

%% Regularization parameters of the fully iterated problem
% These parameters are the starting point of the search below. The optimal
% parameter at a finite budget is typically larger.

R = load(fullfile(toyDir,'RegParameters.mat'));  % bw, Reg
[have,where] = ismember(bw,R.bw);
assert(all(have),'no regularization parameter for the bandwidths %s',mat2str(bw(~have)));
regRef = R.Reg(where);
regRef = regRef(:);

%% LSQR at every time budget
% This is the expensive part of the script. The search for the
% regularization parameter brackets |warmDecades| decades on either side of
% its starting point and bisects |nRefine| times, i.e., it costs
% 3 + 2*nRefine LSQR runs per entry.

nRefine = 2;
warmDecades = 1;

% predicted cost from the time model
tEntry = t0net + slope .* iterBudget;
tPredict = (3 + 2*nRefine) * sum(tEntry(isfinite(tEntry)));
fprintf('\nbudgets from %i to %i iterations, predicted runtime %.1f h\n\n', ...
  min(iterBudget(:),[],'omitnan'),max(iterBudget(:),[],'omitnan'),tPredict/3600);

Reg = nan(nB,nF);        % regularization parameter
ErrorLSQR = nan(nB,nF);  % relative L2-error
tAchieved = nan(nB,nF);  % runtime of the LSQR run
nSolves = zeros(nB,nF);

regPrevBw = nan(1,nF);   % parameters of the previous bandwidth
tSweep = tic;

ws = warning('off','all');
for n = 1:nB
  for k = 1:nF
    it = iterBudget(n,k);

    if ~isfinite(it)
      fprintf('bandwidth %3i, %3gx: no time model, skipped\n',bw(n),timeFactors(k));
      continue
    end

    if it < 1
      % not a single iteration fits into the budget
      [f0,tAchieved(n,k)] = lsqrFit(nodes,values,bw(n),regRef(n),0,lsqrOpts);
      ErrorLSQR(n,k) = relError(f0,fun);
      Reg(n,k) = regRef(n);
      nSolves(n,k) = 1;
      fprintf('bandwidth %3i, %3gx: no iteration, error %.4g\n', ...
        bw(n),timeFactors(k),ErrorLSQR(n,k));
      continue
    end

    % starting point of the search
    if k > 1 && isfinite(Reg(n,k-1))
      reg0 = Reg(n,k-1);
    elseif isfinite(regPrevBw(k))
      reg0 = regPrevBw(k);
    else
      reg0 = regRef(n);
    end

    [Reg(n,k),ErrorLSQR(n,k),tAchieved(n,k),nSolves(n,k)] = ...
      findBestRegularization(nodes,values,fun,bw(n),it,lsqrOpts, ...
      reg0,warmDecades,nRefine);

    regPrevBw(k) = Reg(n,k);

    fprintf(['bandwidth %3i, %3gx: %6i iterations, reg = %.3g, ', ...
      'error %.4g (%.1f of %.1f s)\n'],bw(n),timeFactors(k),it,Reg(n,k), ...
      ErrorLSQR(n,k),tAchieved(n,k) - tVormed,timeFactors(k) * tHAMLS(n));
  end

  % save after every bandwidth, such that an interrupted run keeps its results
  save(fullfile(outDir,'EqualTimeErrorTimeFactors'),'bw','mlsDegree', ...
    'timeFactors','nRefine','warmDecades','Reg','ErrorLSQR','ErrorMLS', ...
    'ErrorProj','iterBudget','tAchieved','tHAMLS','t0net','slope','nSolves');
end
warning(ws);

fprintf('\nfinished after %.2f h and %i LSQR runs\n',toc(tSweep)/3600,sum(nSolves(:)));

%% Write the tables
% The LSQR error jumps between neighbouring bandwidths, because the cost of
% an iteration does: on the machine of the paper an iteration is up to twice
% as expensive for $L \equiv 4 \pmod 8$ as for the neighbouring bandwidths,
% which therefore get fewer iterations out of the same budget. Figure 5 shows
% the raw columns |LSQR_x<f>| at the odd bandwidths only. The table also
% contains the upper envelope |LSQR_x<f>_env| through the local maxima. Which
% bandwidths are expensive depends on the machine, see |slope| in
% |EqualTimeErrorTimeFactors.mat|.

ErrorLSQRenv = nan(nB,nF);
for k = 1:nF
  ErrorLSQRenv(:,k) = maxEnvelope(ErrorLSQR(:,k));
end

errHeader = cell(1,2*nF);
errColumns = nan(nB,2*nF);
for k = 1:nF
  errHeader{2*k-1} = sprintf('LSQR_x%g',timeFactors(k));
  errHeader{2*k} = sprintf('LSQR_x%g_env',timeFactors(k));
  errColumns(:,2*k-1) = ErrorLSQR(:,k);
  errColumns(:,2*k) = ErrorLSQRenv(:,k);
end

itHeader = arrayfun(@(f) sprintf('iter_x%g',f),timeFactors,'UniformOutput',false);

writeTable(fullfile(outDir,'LSQRErrorTimeFactors.txt'), ...
  [{'bandwidth','Projection',sprintf('MLS_Harm_deg%i',mlsDegree)},errHeader], ...
  [bw', ErrorProj, ErrorMLS, errColumns]);

writeTable(fullfile(outDir,'LSQRIterationBudgets.txt'), ...
  [{'bandwidth',sprintf('tHAMLS_deg%i',mlsDegree)},itHeader], ...
  [bw', tHAMLS, iterBudget]);

%% Summary
% We also report how well the time model predicted the actual runtimes. The
% model is calibrated at 0 and 100 iterations but used for up to tens of
% thousands of iterations.

nLast = nB;
kOne = find(timeFactors == 1,1);
kMax = nF;

% the smallest factor at which LSQR reaches the HAMLS error
matchFactor = nan(nB,1);
for n = 1:nB
  k = find(ErrorLSQR(n,:) <= ErrorMLS(n),1);
  if ~isempty(k), matchFactor(n) = timeFactors(k); end
end
nNeverMatched = sum(isnan(matchFactor));

% relative deviation of the actual runtimes from the budgets
tTarget = tHAMLS * timeFactors;
tDev = (tAchieved - tVormed - tTarget) ./ tTarget;
tDev = tDev(isfinite(tDev) & iterBudget >= 1);

fprintf('\nSection 3.2, Figure 5\n');
fprintf('Voronoi areas           %.3f s (subtracted)\n',tVormed);
fprintf('L = %i, HAMLS          %.2f s, error %.3e\n', ...
  bw(nLast),tHAMLS(nLast),ErrorMLS(nLast));
fprintf('L = %i, LSQR at %3gx   %6i iterations, error %.3e\n', ...
  bw(nLast),timeFactors(kOne),iterBudget(nLast,kOne),ErrorLSQR(nLast,kOne));
fprintf('L = %i, LSQR at %3gx   %6i iterations, error %.3e\n', ...
  bw(nLast),timeFactors(kMax),iterBudget(nLast,kMax),ErrorLSQR(nLast,kMax));
fprintf('LSQR never reaches HAMLS up to %gx at %i of %i bandwidths\n', ...
  timeFactors(kMax),nNeverMatched,nB);
if nNeverMatched < nB
  fprintf('median factor at which LSQR reaches HAMLS: %g\n', ...
    median(matchFactor,'omitnan'));
end
if ~isempty(tDev)
  [~,iWorst] = max(abs(tDev));
  fprintf('deviation from the budgets: median %+.0f %%, worst %+.0f %%\n', ...
    100*median(tDev),100*tDev(iWorst));
  if abs(median(tDev)) > 0.1
    warning('The LSQR runs miss their budgets by %+.0f %% in the median.', ...
      100*median(tDev));
  end
end

paperNumbers = struct('nodes',numel(nodes),'degree',mlsDegree, ...
  'timeFactors',timeFactors,'hamlsTimeLast',tHAMLS(nLast), ...
  'hamlsErrorLast',ErrorMLS(nLast),'lsqrIterLast',iterBudget(nLast,:), ...
  'lsqrErrorLast',ErrorLSQR(nLast,:), ...
  'gapAtMaxFactor',ErrorLSQR(nLast,kMax) / ErrorMLS(nLast), ...
  'nNeverMatched',nNeverMatched,'budgetDeviation',median(tDev));

%% A quick look at the results

figure(5); clf
semilogy(bw,ErrorProj,'-','Color',[0.45 0.45 0.45],'LineWidth',1.5); hold on
semilogy(bw,ErrorMLS,'k--','LineWidth',2);
names = {'L2-projection',sprintf('HAMLS (degree %i)',mlsDegree)};
for k = 1:nF
  % from red (short budget) to blue (long budget)
  c = [1 0 0] * (1 - (k-1)/max(nF-1,1)) + [0 0 1] * (k-1)/max(nF-1,1);
  semilogy(bw,ErrorLSQR(:,k),'-','Color',c,'LineWidth',1);
  names = [names, {sprintf('LSQR at %gx',timeFactors(k))}]; %#ok<AGROW>
end
hold off
xlabel('bandwidth L'); ylabel('relative L2-error'); grid on
legend(names,'Location','eastoutside');

%% Helper functions

function [f,t] = lsqrFit(nodes,values,L,alpha,maxit,opts)
% LSQR with maxit iterations and its runtime

tStart = tic;
f = S2FunHarmonic.interpolate(nodes,values,'bandwidth',L, ...
  'regularization',alpha,'maxit',maxit,opts{:});
t = toc(tStart);

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


function env = maxEnvelope(y)
% upper envelope of y through its local maxima
%
% The envelope interpolates linearly between the local maxima and the two
% end points and is never below y. NaN entries remain NaN.

env = nan(size(y));
idx = find(isfinite(y));
if numel(idx) < 3
  env(idx) = y(idx);
  return
end

v = y(idx);
pk = false(size(v));
pk(2:end-1) = v(2:end-1) >= v(1:end-2) & v(2:end-1) >= v(3:end);
pk([1 end]) = true;

e = interp1(idx(pk),v(pk),idx,'linear');
env(idx) = max(e,v);

end


function [reg,err,tBest,nSolves] = findBestRegularization( ...
  nodes,values,fun,L,maxit,opts,reg0,warmDecades,nRefine)
% regularization parameter that minimizes the L2-error for a fixed budget
%
% Syntax
%   [reg,err,tBest,nSolves] = findBestRegularization(nodes,values,fun,L, ...
%     maxit,opts,reg0,warmDecades,nRefine)
%
% Output
%  reg     - best regularization parameter
%  err     - relative L2-error at reg
%  tBest   - runtime of the LSQR run at reg
%  nSolves - number of LSQR runs
%
% The bracket reaches warmDecades decades on either side of reg0 and is
% bisected nRefine times on the logarithmic scale.

nSolves = 0;
best = struct('reg',NaN,'err',Inf,'t',NaN);

m = log10(reg0);
lo = m - warmDecades;
hi = m + warmDecades;

Em = errorFor(10^m);
Elo = errorFor(10^lo);
Ehi = errorFor(10^hi);

% reg0 is the optimum of a neighbouring entry; if the minimum lies outside
% the bracket, we shift the bracket once
if Elo < Em && Elo <= Ehi
  hi = m; m = lo; lo = lo - warmDecades; Em = Elo;
elseif Ehi < Em
  lo = m; m = hi; hi = hi + warmDecades; Em = Ehi;
end

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

reg = best.reg;
err = best.err;
tBest = best.t;

  function e = errorFor(alpha)
    [f,t] = lsqrFit(nodes,values,L,alpha,maxit,opts);
    e = relError(f,fun);
    nSolves = nSolves + 1;
    if e < best.err
      best = struct('reg',alpha,'err',e,'t',t);
    end
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
