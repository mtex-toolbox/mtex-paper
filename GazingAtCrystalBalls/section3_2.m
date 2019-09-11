% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%

%% Section 3.2 Quadrature

clear; close all; home;
tic

%set up cooordinate conventions
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

%% load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% express the master by a harmonic function
cs = loadCIF(cifname);
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);

%% set up the detector

det = detector(400,300,0.56,[0.49,0.45]);

%% simulate a pattern -> Figure 2 a

% generate a random orientation
try
  load ori 
catch
  ori = orientation.rand(cs);
  
  % we would like to have always the same random orientation
  save('ori','ori') 
end

% arguments are flux and background
pattern = det.simulatePattern(master,ori,1000,10);

% and draw it
clf
imagesc(pattern);
axis off
colormap gray
%saveFigure('../pic/detRaw.png')

%% quadrature based approach -> Figure 2 b

tic
pHarm = det.pattern2Fun(pattern,'bandwidth',256,'quadrature','delta',0);
toc

figure(1)
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap gray
%saveFigure('../pic/detQuadrature.png')


%% approximation based approach -> Figure 2 c

tic
pHarm = det.pattern2Fun(pattern,'bandwidth',256,'delta',0,'maxit',5)
toc

figure(1)
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
caxis([-1,100])
colormap gray

%saveFigure('../pic/detApproximation.png')

%% approximation based approach with smooth cutoff -> Figure 2 d

tic
pHarm = det.pattern2Fun(pattern,'bandwidth',256,'delta',0.2,'maxit',5)
toc

figure(1)
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap gray

%saveFigure('../pic/detApproximationCutOff.png')
