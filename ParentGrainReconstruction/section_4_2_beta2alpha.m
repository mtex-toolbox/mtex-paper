% MTEX script for the paper
%
%% Parent grain reconstruction from partially or fully transformed microstructures in MTEX
%
% by Frank Niessen, Tuomo Nyyssönen, Azdiar A. Gazder, Ralf Hielscher:
%
% * to run the script at least MTEX 5.7 has to be installed
% 
%% Section 4.2 - beta to alpha Titanium
%
%% Data import and grain reconstruction

% load the EBSD data
ebsd = mtexdata('alphaBetaTitanium');
ebsd('Ti (alpha)').CS.mineral = 'Alpha';
ebsd('Ti (beta)').CS.mineral = 'Beta';

% grains are calculated with a 1.5° threshold
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',1.5*degree,...
  'removeQuadruplePoints');

% Figure 11 (a) - initial data plot
plot(ebsd('alpha'),ebsd('alpha').orientations)

%% initialize the parent grain reconstructor

job = parentGrainReconstructor(ebsd,grains);

% set orientation relationship to Burgers
job.p2c = orientation.Burgers(job.csParent, job.csChild)

% Figure 11(b) - alpha / alpha boundary misorientation axes
plotIPDF_gB_misfit(job); % requires OR-Tools

%% parent grain reconstruction from triple points

% compute votes
job.calcTPVotes('minFit',2.5*degree,'maxFit',5*degree);

% calcualte parent orientations based on this vote
job.calcParentFromVote('minProb',0.7)

% Figure  11(c) reconstructed parent orientations
plot(job.parentGrains, job.parentGrains.meanOrientation,'linewidth',1.5);

%% grow parent grains

for k = 1:3
  % compute votes
  job.calcGBVotes('p2c','threshold',k * 2.5*degree);
  
  % grow parent grains
  job.calcParentFromVote
end

%% clean reconstruction

% merge grains with similar orientation
job.mergeSimilar('threshold',5*degree);

% merge small inclusions
job.mergeInclusions('maxSize',50);

% Figure 11(d) plot reconstructed EBSD
plot(job.ebsd('beta'),job.ebsd('beta').orientations)

% plot reconstrcuted grain boundaries on top
hold on
plot(job.grains.boundary,'lineWidth',3)
hold off

