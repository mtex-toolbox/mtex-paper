%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Optimal Sampling on SO(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run with MatlabR2024b

clear

cs = crystalSymmetry.load('Mg-Magnesium.cif');
cs = cs.properGroup;
odf = fibreODF(cs.cAxis, vector3d.Z,'halfwidth',20*degree);
odf = SO3FunHarmonic(odf);

% Sampling Sizes
M = [8,16,32,64,128,256,512,1024,2048,4096,8192];

% harmonic degree up to which the sample is optimized, i.e. the orientation
% resolution the sample is required to reproduce. It has to be stated
% explicitly, since the discrepancy below is only meaningful with respect to
% the same bandwidth. 32 is the default of optimalSample.
bw = 32;


%% Compute Optimal Sampling

maxIter = 10000;
tol = 0.0001*degree;

% compute optimal sets of size M
% We start by a discrete sample, since we can not compute equispaced grids of specific length
oriOpt = {};
for i = 1:length(M)
  fprintf(['Number of Points: ',num2str(M(i)),'\n'])
  rng(0)
  oriOpt{i} = odf.discreteSample(M(i),'compact','bandwidth',bw,'maxIter',maxIter,'tol',tol,'steepestDescent');
  save('OptimalOri','oriOpt','M','bw')
end


%% Compute Random Sampling

rng(0)
oriRand = {};
for i = 1:length(M)
  for j=1:100
    fprintf(['Number of Points: ',num2str(M(i)),'  Iteration: ',num2str(j),'\n'])
    oriRand{i,j} = odf.discreteSample(M(i));
  end
  save('RandomOri','oriRand','M')
end



%% Compute Error

% Every sample is measured in two ways.
%
% (1) The relative L2 error of an ODF reconstructed from the sample by kernel
% density estimation, with the halfwidth tuned individually for every sample
% so that no sample is penalized by a badly chosen kernel. This is the error
% a user of the sample sees, but it measures the reconstruction and not the
% sample alone, and it is dominated by the smoothing bias of the kernel.
%
% (2) The kernel discrepancy of the sample itself, see discrepancySO3, i.e.
% exactly the functional E that the optimal sampling minimizes. It needs no
% smoothing parameter and it bounds the quadrature error of the orientation
% set for every property of bandwidth bw.

load('RandomOri')
% Kernel Density Estimation and Error for Random Sampling
for i=1:length(M)
  for j=1:100
    fprintf(['Number of Points: ',num2str(M(i)),'  Iteration: ',num2str(j),'\n'])
    ori = oriRand{i,j}; ori.CS = cs; % TODO: Bug fix resymmetrise
    errorFun = @(hw) norm(odf-calcDensity(ori,'halfwidth',hw*degree))/norm(odf);
    [hwOpt, errOpt] = fminbnd(errorFun, 5, 25,optimset('TolX',0.1));
    hw_RandSampl(i,j) = hwOpt;
    Erel_RandSampl(i,j) = errOpt;
    Ediscr_RandSampl(i,j) = discrepancySO3(odf,ori,[],'bandwidth',bw);
  end
  save('RandomSamplingError','Erel_RandSampl','hw_RandSampl','Ediscr_RandSampl','M','bw')
end


% Kernel Density Estimation and Error for Optimal Sampling
setMTEXpref('maxSO3Bandwidth',92)
load('OptimalOri')
for i=1:length(M)
  fprintf(['Number of Points: ',num2str(M(i)),'\n'])
  ori = oriOpt{i}; ori.CS = cs; % TODO: Bug fix resymmetrise
  errorFun = @(hw) norm(odf-calcDensity(ori,'halfwidth',hw*degree))/norm(odf);
  [hwOpt, errOpt] = fminbnd(errorFun, 2.5, 25,optimset('TolX',0.1));
  hw_OptSampl(i) = hwOpt;
  Erel_OptSampl(i) = errOpt;
  Ediscr_OptSampl(i) = discrepancySO3(odf,ori,[],'bandwidth',bw);
  save('OptimalSamplingError','Erel_OptSampl','hw_OptSampl','Ediscr_OptSampl','M','bw')
end



%% Plotting the Data

% 3d plot of odf
figure(1)
plot3d(odf,'AxisAngle')
h = gcf(); h.Children(2).Visible = 'off';
set(h,'position',[10,10,1000,1000])
drawnow
exportgraphics(h,'ODF3d.png','Resolution',300);
drawnow
close(h)

i0 = 8;

load('OptimalOri')
load('OptimalSamplingError')
fprintf(['Plot Density Estimation and Optimal Sampling for M=',num2str(M(i0)),'orientations.\n'])

% 3d plot of optimal points
figure(1)
scatter(oriOpt{i0},'AxisAngle','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','b')
h = gcf(); h.Children(2).Visible = 'off';
set(h,'position',[10,10,1000,1000])
drawnow
exportgraphics(h,'Compactification3d.png','Resolution',300);
drawnow
close(h)

% 3d plot of Density Estimation
figure(1)
Rec_odf = calcDensity(oriOpt{i0},'halfwidth',hw_OptSampl(i0)*degree);
plot3d(Rec_odf,'AxisAngle')
h = gcf(); h.Children(2).Visible = 'off';
set(h,'position',[10,10,1000,1000])
drawnow
exportgraphics(h,'DensityEstimation3d.png','Resolution',300);
drawnow
close(h)





%% Plotting the Numerics

% load data
load('RandomSamplingError.mat')
load('OptimalSamplingError')

figure(1)
% Boxplots
boxchart( repmat((1:11)', size(Erel_RandSampl,2), 1), Erel_RandSampl(:), 'BoxFaceAlpha', 0.25);
hold on
% Error plots
plot(1:11, Erel_OptSampl, '-o','LineWidth', 2);
plot(1:11, mean(Erel_RandSampl,2), '-o','LineWidth', 2);
hold off

% Plot Properties
xticks(1:11)
xticklabels(string(M))
xlabel('Number of Points')
ylabel('Relative error')
legend('Random sampling','Optimal sampling','Random sampling mean','Location', 'best')
grid on



%% Precompute and save the data for tikz

load('RandomSamplingError.mat')
load('OptimalSamplingError')

% Data for line segments
ErelOpt = Erel_OptSampl(:);
MeanRand = mean(Erel_RandSampl, 2);

DiscrOpt = Ediscr_OptSampl(:);
DiscrMeanRand = mean(Ediscr_RandSampl, 2);

% Data for Boxplots, for both error measures
for k = 1:length(M)
    y = Erel_RandSampl(k,:);
    q = quantile(y, [0.25, 0.5, 0.75]);

    Q1(k,1) = q(1);
    Median(k,1) = q(2);
    Q3(k,1) = q(3);
    lowerBound = q(1) - 1.5*(q(3)-q(1));
    upperBound = q(3) + 1.5*(q(3)-q(1));
    inliers = y >= lowerBound & y <= upperBound;
    LowerWhisker(k,1) = min(y(inliers));
    UpperWhisker(k,1) = max(y(inliers));

    y = Ediscr_RandSampl(k,:);
    q = quantile(y, [0.25, 0.5, 0.75]);

    DiscrQ1(k,1) = q(1);
    DiscrMedian(k,1) = q(2);
    DiscrQ3(k,1) = q(3);
    lowerBound = q(1) - 1.5*(q(3)-q(1));
    upperBound = q(3) + 1.5*(q(3)-q(1));
    inliers = y >= lowerBound & y <= upperBound;
    DiscrLowerWhisker(k,1) = min(y(inliers));
    DiscrUpperWhisker(k,1) = max(y(inliers));
end

% Export
% The column names are the contract with Kapitel/OptimalSampling.tex, which
% reads this file by \pgfplotstableread.
NumPoints = M';
T = table( NumPoints, ErelOpt, MeanRand, Q1, Median, Q3, LowerWhisker, UpperWhisker, ...
  DiscrOpt, DiscrMeanRand, DiscrQ1, DiscrMedian, DiscrQ3, ...
  DiscrLowerWhisker, DiscrUpperWhisker);
writetable(T, 'SamplingError.txt', 'Delimiter', 'tab');



