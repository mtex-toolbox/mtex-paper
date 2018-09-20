% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 2.7 Morphological Filters and Inpaiting
%

% where to save the image files
pname = ['..' filesep 'pic' filesep 'sim' filesep];
plotx2east

% de la Vallee Poussin distributed noise
ebsd = simEBSD('poussin',2*degree,'saltPepper',0.05);

%% outliers identified as one pixel grains - Fig. 2.8a

[grains,ebsd.grainId] = calcGrains(ebsd)

plot(ebsd,angle(ebsd.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
caxis([0,18])

hold on
plot(grains.boundary,'lineWidth',2)
hold off


%% one pixel grains removed - Fig. 2.9b

ebsd(grains(grains.grainSize <= 1)) = [];

plot(ebsd,angle(ebsd.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
caxis([0,18])

%saveFigure('../pic/sim/grainsRemoved.png')

%% mean filter - Fig. 2.8c

F = meanFilter; F.weights = ones(5);

ebsd_smoothed = smooth(ebsd,F,'fill')
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure('../pic/sim/meanGrain.png')

%% median filter - Fig. 2.8d

F = medianFilter; F.numNeighbours = 3;

ebsd_smoothed = smooth(ebsd,F,'fill')
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure('../pic/sim/medianGrain.png')

%% Kuwahara filter - Fig. 2.8e

F = KuwaharaFilter; F.numNeighbours = 3;

ebsd_smoothed = smooth(ebsd,F,'fill')
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')
%saveFigure('../pic/sim/kuwaharaGrain.png')

%% smoothing spline - Fig. 2.8e

F = splineFilter; F.alpha = 5;

ebsd_smoothed = smooth(ebsd,F,'fill')
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure('../pic/sim/splineGrain.png')

%% total variation - Fig. 2.8f

F = halfQuadraticFilter; F.alpha = 0.5;

ebsd_smoothed = smooth(ebsd,F,'fill')
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure('../pic/sim/halfquadGrain.png')
