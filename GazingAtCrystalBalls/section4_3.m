% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%

%% Section 4.3 Peak Detection

clear; close all; home;
tic

%% Load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));

% as hamonic expansion
cs = loadCIF(cifname);
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);

%% define a band profile

% modified Gaussian profile
profile = @(x) 1.5*exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-87*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-93*degree).^2./(2*degree).^2);

% expand into a Legendre series
profileHarm = S2Kernel.quadrature(profile);

%% simulate a noisy pattern

load ori

det = detector(400,300,0.56,[0.49,0.45]);

% define a cut off function for the detector
mask = det.S2CutOffMask(0);
maskHarm = S2FunHarmonic.quadrature(mask);
correction = max(maskHarm.radon,0.1);

pattern = det.simulatePattern(master,ori,100,10);

% plotting
imagesc(pattern);
colormap gray; axis off
% saveFigure('../pic/noisyMap.png')

%%  approximate the pattern by a harmonic function

% approximate pattern by a harmonic function
pHarm = det.pattern2Fun(pattern,'bandwidth',256,'quadrature','delta',0.1);

plot(pHarm,'upper')
colormap gray

%% compute corrected spherical convolution 

% 
RadonPHarm = conv(pHarm,profileHarm)./ correction;
RadonPHarm.bandwidth = 64;

plotx2east
figure(1)
plot(RadonPHarm,'pcolor','resolution',0.25*degree,'upper')
colormap gray


%% Peak detection -> Figure 5b

peaks = S2PeakDetection(RadonPHarm,det,15,-0.11);

% mark found normals
annotate(peaks,'markerSize',10,'antipodal','MarkerEdgeColor','blue','MarkerFaceColor','none')

% mark theoretical band normals
h = Miller({1,0,0},{1,1,0},{1,1,1},{2,1,0},{3,1,0},{2,1,1},{3,2,1},{3,3,2},cs);
annotate(ori*h([1,2,6]).symmetrise,'markerSize',15,'antipodal','MarkerEdgeColor','r','MarkerFaceColor','none','Marker','o')
%saveFigure('../pic/peakDetection.png')

%% Annotate bands in sperical Kikuchi pattern -> Figure 5a

figure(2)

% plotting
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap gray

circle(peaks,'linecolor',[ 0.3 0.3 1])

%saveFigure('../pic/bandsSphere.png')
