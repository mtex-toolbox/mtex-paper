% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 3.2 Gradient and Curvature Tensor Estimation
%


%% The raw data in ipf colors - Fig. 3.4a

% some plotting conventions
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');


% import the h5 file to MTEX
% header is a struct which contains all the microscope data
ebsd = loadEBSD_h5('4_HREBSD_1degreeNP.h5');


% grain reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',4*degree);

% remove small grains - 
grains = grains(grains.area >= 10); %remove grains

% and throw away these measurements from the ebsd data set
ebsd = ebsd(grains);

% redo grain reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',4*degree);

% and smooth the grain boundaries a bit to avoid staircasing effect
% this should be done before denoising the ebsd data
grains = smooth(grains,3);


plot(ebsd('fe'),ebsd('fe').orientations,'micronbar','off')
hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('ebsdRAW.png')


%% raw data in axis angle colors - Fig. 3.4b

ipfKey = axisAngleColorKey(ebsd);
ipfKey.oriRef = grains(ebsd('fe').grainId).meanOrientation;

plot(ebsd('fe'),ipfKey.orientation2color(ebsd('fe').orientations),'micronbar','off')
hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdAxisAngle.png')

%% data denoising - Fig. 3.4c

F = halfQuadraticFilter;
F.alpha = 0.1;

ebsd_smoothed = smooth(ebsd,F,'fill',grains); 

plot(ebsd_smoothed('fe'),ebsd_smoothed('fe').orientations)
hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdSRAW.png')

%% denoised data in axis angle colors - Fig. 3.4d

ipfKey = axisAngleColorKey(ebsd);
ipfKey.oriRef = grains(ebsd_smoothed('fe').grainId).meanOrientation;

plot(ebsd_smoothed('fe'),ipfKey.orientation2color(ebsd_smoothed('fe').orientations),'micronbar','off')
hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdSAxisAngle.png')

%% lattice curvature tensor of the noisy data - Fig. 3.5a 

ebsd = ebsd('fe').gridify;

kappa = ebsd.curvature;

plot(ebsd,abs(kappa{1}),'micronbar','off','outermargin',20)
set(gca,'ColorScale','log')
caxis([0.001,0.1])
mtexColorMap LaboTeX
mtexColorbar('fontSize',30)

hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdKappa11.png')

%% lattice curvature tensor of the denoised data - Fig. 3.5b

ebsdS = ebsd_smoothed('indexed').gridify;

kappa = ebsdS.curvature;

plot(ebsdS,abs(kappa{1}),'micronbar','off')
set(gca,'ColorScale','log')
caxis([0.001,0.1])
mtexColorMap LaboTeX
mtexColorbar('fontSize',30)

hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdSKappa11.png')

