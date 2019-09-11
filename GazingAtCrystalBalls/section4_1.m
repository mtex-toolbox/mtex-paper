% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%

%% Section 4.1 The Spherical Radon Transform

clear; close all; home;
setMTEXpref('maxBandwidth',256)

%% Load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));

% as hamonic expansion
cs = loadCIF(cifname);
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);


%% Radon transform of the master pattern -> Figure 3a

plot(masterHarm.radon,'pcolor','resolution',0.25*degree,'upper','complete')
colormap gray

% highlight some plane normals
h = Miller({1,0,0},{1,1,0},{1,1,1},{2,1,1},{3,1,0},cs)

colors = vega10(5);

hold on
for i = 1:length(h)
  plot(h(i).symmetrise,'MarkerSize',25,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',26,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',27,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',23,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',24,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
end
hold off

%saveFigure('../pic/radonMaster.png')

