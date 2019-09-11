s% Script file for the publication
%
%% Gazing at crystal balls 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%
%% Section 3.1 Quadrature

clear; close all; home;
tic

%% Load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));

%% spherical projection of the master pattern -> Figure 1 a

plot(master,'upper','pcolor','resolution',0.25*degree,'complete')
colormap gray

%saveFigure('../pic/master.png')

%% Quadrature based approximation -> Figure 1 b,c,d

% as hamonic expansion
cs = loadCIF(cifname);
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',256);

plot(masterHarm,'upper','pcolor','resolution',0.25*degree,'complete')
colormap gray

%saveFigure('../pic/master256.png')