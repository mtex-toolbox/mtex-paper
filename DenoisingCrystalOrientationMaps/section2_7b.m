% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 2.7 Morphological Filters and Inpaiting
%

close all; plotx2east
mtexdata forsterite

%restrict to a subregion
ebsd = ebsd(inpolygon(ebsd,[10 4 5 3]*10^3));


%% the raw geological data - Fig. 2.9a

ipfKey = ipfHSVKey(ebsd('Fo')); ipfKey.inversePoleFigureDirection = vector3d(1,1,0);
plot(ebsd('Fo'),ipfKey.orientation2color(ebsd('Fo').orientations),'micronbar','off')
hold on
plot(ebsd('En'),ebsd('En').orientations)
plot(ebsd('Di'),ebsd('Di').orientations)
hold off

%saveFigure('../pic/ebsdGrainRaw.png')

%% traditional grain reconstruction - Fig. 2.9b

[grains,ebsd.grainId] = calcGrains(ebsd,'angle',10*degree)

% start overide mode
hold on

% plot the boundary of all grains
plot(grains.boundary,'linewidth',2)

% stop overide mode
hold off

%saveFigure('../pic/ebsdGrain1.png')

%% advanced grain reconstruction - Fig. 2.9c 

plot(ebsd('Fo'),ipfKey.orientation2color(ebsd('Fo').orientations),'micronbar','off')
hold on
plot(ebsd('En'),ebsd('En').orientations)
plot(ebsd('Di'),ebsd('Di').orientations)

[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',10*degree)

ebsd(grains(grains.grainSize<3)) = [];

[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',10*degree)

grains = smooth(grains,5);

% start overide mode

% plot the boundary of all grains
plot(grains.boundary,'linewidth',2)

% stop overide mode
hold off

%saveFigure('../pic/ebsdGrain2.png')


%% axis / angle to grain meanorientation map - Fig. 2.10a
% For comparison
ipfKey = axisAngleColorKey(ebsd('Fo'));
ipfKey.oriRef = grains(ebsd('fo').grainId).meanOrientation;
ipfKey.maxAngle = 2.5*degree;

color = ipfKey.orientation2color(ebsd('Fo').orientations);
plot(ebsd('Fo'),color,'micronbar','off')

hold on
ipfKey.oriRef = grains(ebsd('En').grainId).meanOrientation;

plot(ebsd('En'),ipfKey.orientation2color(ebsd('En').orientations))

% plot boundary
plot(grains.boundary.reorder,'linewidth',3)
plot(grains('En').boundary.reorder,'lineWidth',3)
hold off

%saveFigure('../pic/FoAxisAngle.png')


%% smooth the data 

% select the filter of your choice
%F = meanFilter; F.weights = ones(5);
%F = medianFilter; F.numNeighbours = 3;
%F = splineFilter
%F = KuwaharaFilter; F.numNeighbours = 3;
%F = medianFilter; F.numNeighbours = 1;
F = halfQuadraticFilter; F.alpha = 0.01;


ebsd_smoothed = smooth(ebsd('indexed'),F,'fill',grains)

%% denoised orientation maps - Fig. 2.10 (b) - (f)

ipfKey = axisAngleColorKey(ebsd_smoothed('Fo'));
ipfKey.oriRef = grains(ebsd_smoothed('fo').grainId).meanOrientation;
ipfKey.maxAngle = 2.5*degree;

color = ipfKey.orientation2color(ebsd_smoothed('Fo').orientations);
plot(ebsd_smoothed('Fo'),color,'micronbar','off')

hold on
ipfKey.oriRef = grains(ebsd_smoothed('En').grainId).meanOrientation;

plot(ebsd_smoothed('En'),ipfKey.orientation2color(ebsd_smoothed('En').orientations),'figSize','tiny')

% plot boundary
plot(grains.boundary.reorder,'linewidth',3)

hold off

%saveFigure('../pic/FoHalfQuad.png')
%saveFigure('../pic/FoSpline.png')
%saveFigure('../pic/FoMean.png')
%saveFigure('../pic/FoMedian.png')
%saveFigure('../pic/FoKuwahara.png')
