% MTEX script for the paper
%
%% Parent grain reconstruction from partially or fully transformed microstructures in MTEX
%
% by Frank Niessen, Tuomo Nyyssönen, Azdiar A. Gazder, Ralf Hielscher:
%
% * to run the script at least MTEX 5.7 has to be installed
% 
%% Section 4.1 - gamma to alpha in lath martensite
% 
%
%% Data Import and Grain Reconstruction

% import the data and
% rename Iron bcc (old) to 'AlphaP' and 'Iron fcc' to 'Gamma'
ebsd = mtexdata('martensite');
ebsd('Iron bcc (old)').CS.mineral = 'AlphaP';
ebsd('Iron fcc').CS.mineral = 'Gamma';

% calculate grains with a 3° threshold
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 3*degree);

% remove EBSD data corresponding to one pixel grains
ebsd(grains(grains.grainSize < 3)) = [];

% recalculate the grains from the remaining data ...
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',3*degree);

% smooth grain boundaries
grains = smooth(grains,5);

%% Figure 1 - the initial map

% set default plotting convention
plotx2east; plotzOutOfPlane;

plot(ebsd,ebsd.orientations)
hold on
plot(grains.boundary)
hold off

%% Figure 2a - disorientation histogram before and after OR refinement

% set up the parent grain reconstructor
job = parentGrainReconstructor(ebsd,grains);

% initial OR guess: Kurdjumow-Sachs
job.p2c = orientation.KurdjumovSachs(job.csParent, job.csChild);

fit = job.calcGBFit;
close all
histogram(fit./degree,'DisplayName','K-S');

% refine OR based on the fit with boundary misorientations
job.calcParent2Child;
hold on
fit = job.calcGBFit;
histogram(fit./degree,'DisplayName','refined');
hold off

%% Figure 2b - variant map

% display variants in the (001) pole figure
job.plotVariantPF(Miller({0,0,1},job.csChild),'reduced')


%% Figure 3 - boundary misfit map

% compute the misfit for all child to child grain neighbours
[fit,c2cPairs] = job.calcGBFit;

% select grain boundary segments by grain ids
[gB,pairId] = job.grains.boundary.selectByGrainId(c2cPairs);

% plot the child phase
plot(ebsd('alphaP'),ebsd('alphaP').orientations,'figSize','large','faceAlpha',0.5)

% and on top of it the boundaries colorized by the misfit
hold on;

% scale fit between 0 and 1 - required for edgeAlpha
plot(gB, 'edgeAlpha', (fit(pairId) ./ degree - 2.5)./2 ,'linewidth',2);
hold off


%% Figure 4 - graph clusters

% build graph and cluster it using the MCL algorithm
job.calcGraph('threshold',2.5*degree,'tolerance',1.5*degree);
job.clusterGraph('inflationPower',1.6);

% visualize the clusters -> requires the ORTools addon
plotMap_clusters(job,vector3d.Z,'linewidth',2);

%% Figure 5 - parent orientations computed from the graph clusters

% calculate the parent orientations
job.calcParentFromGraph;

% plot reconstructed parent microstructure
plot(job.parentGrains, job.parentGrains.meanOrientation, 'linewidth',2);

%% Figure 6 - Misfit between reconstruction 

plot(job.grains, job.grains.fit./degree,'linewidth',2);
setColorRange([0,5]);
colormap('viridis')
mtexColorbar;

%% Figure 7 - revert parent grains that are too small or have bad fit

job.revert(job.grains.fit > 5*degree | job.grains.clusterSize < 10)

% plot the remaining grains
plot(job.parentGrains, job.parentGrains.meanOrientation, 'linewidth',2)

%% Figure 8 - grow parent grains

for k = 1:3 

  % compute votes
  job.calcGBVotes('p2c','threshold',2.5*degree);
  
  % compute parent orientations from votes
  job.calcParentFromVote

end

% plot reconstructed parent microstructure
plot(job.parentGrains,job.parentGrains.meanOrientation,'linewidth',2)

%% Figure 9 - merge and clean up reconstructed grains

% merge grains with similar orientation
job.mergeSimilar('threshold',7.5*degree);

% merge small inclusions
job.mergeInclusions('maxSize',50);

plot(job.parentGrains,job.parentGrains.meanOrientation,'linewidth',2)

%% Figure 10 - parent EBSD map

% compute the parent EBSD
parentEBSD = ebsd.parentEBSD;

% plot the reconstructed EBSD data
plot(parentEBSD, parentEBSD.orientations);

% with the reconstructed beta grain boundaries
hold on; 
plot(job.grains.boundary,'lineWidth',3)
hold off


%% Figure 13 - variant analysis

% compute variantId and packetId for each priorGrain
job.calcVariants;

% plot packet Id
color = ind2color(job.grainsPrior.packetId);
plot(job.grainsPrior,color)

%%

% variant Id
plot(job.grainsPrior,job.grainsPrior.variantId)
mtexColorMap(jet(24))

%% Figure 15: refine gamma twins 

% check the gamma grains interactively - requires ORTools addon
grainClick(job,parentEBSD);

%%

% check for twins - requires ORTools addon
grainClick(job,parentEBSD,'parentTwins');

