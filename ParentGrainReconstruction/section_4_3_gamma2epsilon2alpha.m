% MTEX script for the paper
%
%% Parent grain reconstruction from partially or fully transformed microstructures in MTEX
%
% by Frank Niessen, Tuomo Nyyssönen, Azdiar A. Gazder, Ralf Hielscher:
%
% * to run the script at least MTEX 5.7 has to be installed
% 
%% Section 4.3 - Gamma to Epsilon to Alpha TWIP/TRIP steel
%
%% data import and grain reconstruction

% load data from file and rename phase names
ebsd = EBSD.load('TRWIP_CR10_E7_1C.cpr');
ebsd('Iron fcc').CS.mineral = 'Gamma';
ebsd('Iron bcc').CS.mineral = 'AlphaP';
ebsd('Epsilon').CS.mineral = 'Epsilon';

% reconstruct grains
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 3*degree);

% EBSD data in small grains are removed
ebsd(grains(grains.grainSize < 3)) = [];

% recalculate the grains from the remaining data
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',3*degree);

% smooth the grain boundaries
grains = smooth(grains,5);

%% Figure 12(a) - phase map

plot(grains)

%% Figure 12(c) - gamma orientations

plot(ebsd('gamma'),ebsd('gamma').orientations)

hold on
plot(grains.boundary)
hold off

%% Figure 12(c) - epsilon orientations

plot(ebsd('epsilon'),ebsd('epsilon').orientations)

hold on
plot(grains.boundary)
hold off

%% Figure 12(d) - alpha' orientations

plot(ebsd('alphaP'),ebsd('alphaP').orientations)

hold on
plot(grains.boundary)
hold off

%% Epsilon to AlphaP Reconstruction

% define reconstruction job
job = parentGrainReconstructor(ebsd,grains);

% set OR to Burgers and refine
job.p2c = orientation.Burgers(ebsd('epsilon').CS,ebsd('alphaP').CS);
job.calcParent2Child('p2c')

% grow epsilon phase
for k = 1:3
  job.calcGBVotes('p2c','threshold',k*2.5*degree);
  job.calcParentFromVote;
end
job.mergeSimilar('threshold',7.5*degree);

% Figure 12(f) - reconstructed epsilon orientations
plot(job.ebsd('epsilon'),job.ebsd('epsilon').orientations)
hold on
plot(job.grains.boundary)
hold off

%% Gamma to Epsilon Reconstruction

% define reconstruction job
job2 = parentGrainReconstructor(job.ebsd,job.grains);

% set OR to Shoji Nishiyama and refine
job2.p2c = orientation.ShojiNishiyama(ebsd('gamma').CS,ebsd('epsilon').CS);
job2.calcParent2Child('p2c')

% grow gamma phase
for k = 1:3
  job2.calcGBVotes('p2c','threshold',k*2.5*degree);
  job2.calcParentFromVote;
end
job2.mergeSimilar('threshold',7.5*degree);

% Figure 12(e) - final reconstruction
plot(job2.parentEBSD, job2.parentEBSD.orientations)
hold on
plot(job2.grains.boundary)
hold off
