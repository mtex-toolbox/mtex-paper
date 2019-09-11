% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%

%% Section 4.2 Spherical convolution and band localisation


clear; close all; home;
tic

% in this section we require also the chebfun package
% https://www.chebfun.org/
addpath ~/repo/chebfun/


%% Load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));

% as hamonic expansion
cs = loadCIF(cifname);
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);


%% compute band profiles -> Figure 3 b

% select some plane normals
h = Miller({1,0,0},{1,1,0},{1,1,1},{2,1,1},{3,1,0},cs);

% compute band profiles by averaging the spherical pattern about all
% rotations around the plane normals
for i = 1:length(h)
  [~,psi{i}] = masterHarm.symmetrise(h(i))
  psi{i}.A(1) = 0;
end

figure
for i = 1:length(h)
 
  plot(psi{i},'DisplayName',char(h(i)),'linewidth',3)
  hold on
  
end

xlim([80,100])
hold off
legend

%saveFigure('../pic/profiles.pdf')

%% spherical convolution with profile (211) -> Figure 4 a

% only consider the profile within a range of 5 degree
% and compute its Legendre series expansion
profileHarm211 = S2Kernel.quadrature(@(t) (abs(t)<cos(85*degree)) .* psi{4}.eval(t))

% plot the profile
figure
plot(profileHarm211,'linewidth',2)
set(gca,'FontSize',20)
xlim([80,100])
axis off;
%saveFigure('../pic/profile3.pdf')

% compute the spherical convolution with the (211) profile
rF211 = conv(masterHarm,profileHarm211);

% plot the spherical convolution
figure;
plot(rF211,'resolution',0.25*degree,'pcolor','upper','complete')
mtexColorMap black2white

%saveFigure('../pic/profile3Conv.png')

%% spherical convolution with profile (310) -> Figure 4 b

% only consider the profile within a range of 5 degree
% and compute its Legendre series expansion
profileHarm310 = S2Kernel.quadrature(@(t) (abs(t)<cos(85*degree)) .* psi{5}.eval(t))

% plot the profile
plot(profileHarm310,'linewidth',2)
set(gca,'FontSize',20)
xlim([80,100])
axis off;
%saveFigure('../pic/profile4.pdf')

% compute the spherical convolution with the (310) profile
rF310 = conv(masterHarm,profileHarm310);

% plot the spherical convolution
figure;
plot(rF310,'resolution',0.25*degree,'pcolor','upper','complete')
mtexColorMap black2white

%saveFigure('../pic/profile4Conv.png')

%% spherical convolution with a Gaussian profil -> Figure 4 c

% define a Gaussian profile
profileGauss = @(x) exp(-(acos(x)-90*degree).^2./(2*degree).^2) ;

% expand it into a Legendre series
profileHarmGauss = S2Kernel.quadrature(profileGauss);

% plot the Gauss profile
figure
plot(profileHarmGauss,'linewidth',2)
set(gca,'FontSize',20)
xlim([80,100])
axis off;
%saveFigure('../pic/profile1.pdf')

% compute the spherical convolution with the Gaussian profile
rFGauss = conv(masterHarm,profileHarmGauss);

% plot the spherical convolution
figure;
plot(rFGauss,'resolution',0.25*degree,'pcolor','upper','complete')
caxis([0.0155 0.0175])
mtexColorMap black2white

%saveFigure('../pic/profile1Conv.png')


%% spherical convolution with a modified Gaussian profil -> Figure 4 d

% define a modified Gaussian profile
profileModGauss = @(x) exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-86*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-94*degree).^2./(2*degree).^2)-0.5;

% expand it into a Legendre series
profileHarmModGauss = S2Kernel.quadrature(profileModGauss);

% plot the profile
figure
plot(profileHarmModGauss,'linewidth',2)
set(gca,'FontSize',20)
xlim([80,100])
axis off;
%saveFigure('../pic/profile2.pdf')

% compute the spherical convolution with the Gaussian profile
rFModGauss = conv(masterHarm,profileHarmModGauss);

% plot the spherical convolution
figure;
plot(rFModGauss,'resolution',0.25*degree,'pcolor','upper','complete')
mtexColorMap black2white
%saveFigure('../pic/profile2Conv.png')
