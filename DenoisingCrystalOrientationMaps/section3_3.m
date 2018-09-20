% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 3.3 Dislocation Density Estimation
%
% please run section 3.2 first


%% compute GND from noisy data

% define bcc dislocation systems
dS = dislocationSystem.bcc(ebsd.CS)

% compute GND - this will take a while
gnd = ebsd.calcGND(dS);

%% plot GND - Fig. 3.6a

plot(ebsd,gnd,'micronbar','off')
set(gca,'ColorScale','log')
caxis([1e13,1e15])
mtexColorMap LaboTeX
mtexColorbar('fontSize',30)

hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdGND.png')

%% compute GND from denoised data

% this will take a while
gndS = ebsdS.calcGND(dS);

%%  plot GND - Fig. 3.6b

plot(ebsdS,gndS,'micronbar','off')
set(gca,'ColorScale','log')
caxis([1e13,1e15])
mtexColorMap LaboTeX
mtexColorbar('fontSize',30)

hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdSGND.png')

%%   plot GND data from high resolution EBSD - Fig. 3.6c

load('4_HREBSD_1degreeNP_HREBSD.mat')
plot(ebsd_smoothed,GND.total,'micronbar','off')
set(gca,'ColorScale','log')
caxis([1e13,1e15])
mtexColorMap LaboTeX
mtexColorbar('fontSize',30)

hold on
plot(grains.boundary,'lineWidth',2)
hold off

%saveFigure('../pic/ebsdCGND.png')
