% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%
%% Section 5 Spherical cross correlation based orientation determination


%% load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));

% as hamonic expansion
cs = loadCIF(cifname);
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);

%% set up the detector

det = detector(400,300,0.56,[0.49,0.45]);

% load ori.mat
% ori = orientation.rand(cs);
ori = orientation('Euler',[30 40 50]*degree);

pattern = det.simulatePattern(master,ori,10000,0);

% approximate the pattern by a harmonic function -
pHarm = det.pattern2Fun(pattern,'bandwidth',128,'quadrature');

%% the uncorrected spherical cross correlation function

% the harmonic bandwidth for the cross correlation function
bw = 128;
setMTEXpref('maxBandwidth',bw)

% compute spherical convolution
xcor = conv(masterHarm,pHarm,'bandwidth',bw);

%% one dimensional section of the spherical cross correlation function -> Figure 7a

figure

omega = linspace(-15,15)*degree;
rot = rotation.byAxisAngle(xvector,omega);
plot(omega./degree,xcor.eval(rot * ori))

%ylim([-0.01,0.2])
xlabel('misorientation angle')
ylabel('cross correlation')
legend('N=128','N=64','N=32')

%saveFigure('../pic/xcor1d.pdf')


%% plot as phi1 section -> Figure 7b

%plot(xcor,'sigma',(0:10:80)*degree,'pcolor','resolution',0.5*degree,'sigma','sections',12)
plot(xcor,'phi1',30*degree,'pcolor','resolution',0.25*degree,'figSize','small')
mtexColorMap LaboTeX

%saveFigure('../pic/xcorPhi1Pattern.png')

%% the correction function -> Figure 7c

% define a cut off function for the detector
mask = det.S2CutOffMask(0.05);

% compute the auto correlation of the pattern with the cutoff function
correction = 1./sum(mask) * conv(masterHarm, mask,'bandwidth',bw) ;

plot(correction,'phi1',30*degree,'pcolor','resolution',0.25*degree,'figSize','small')
mtexColorMap LaboTeX
caxis([0.5120 0.53])
%saveFigure('../pic/xcorPhi1Correction.png')

%% the corrected correlation function -> Figure 7d

% correct for pattern shape
xcor = conv(masterHarm,pHarm,'bandwidth',bw);
xcor = FourierODF(xcor - sum(pHarm)*correction);

%plot(xcor,'sigma',(0:10:80)*degree,'pcolor','resolution',0.5*degree)
plot(xcor,'phi1',30*degree,'pcolor','resolution',0.25*degree,'figSize','small')

mtexColorMap LaboTeX

%saveFigure('../pic/xcorPhi1PatternCorrected.png')

%% peak detection with local refinement - the parameters

% the harmonic bandwidth for the cross correlation function
bw = 64;
setMTEXpref('maxBandwidth',bw)

% the global grid for peak detection
S3G = equispacedSO3Grid(cs,'resolution',2.5*degree);

% the local grid for peak detection
localGrid = localOrientationGrid(specimenSymmetry,specimenSymmetry,...
  2.5*degree,'resolution',0.2*degree);

%% set up the detector and the correction function

det = detector(400,300, 0.57,[0.5,0.45]);

% define a cut off function for the detector
mask = det.S2CutOffMask(0.05);

% compute the auto correlation of the pattern with the cutoff function
correction = 1/sum(mask) * conv(masterHarm, mask,'bandwidth',bw) ;

%% Loop for multiple runs
% set maxit=1 for a single run

tic
maxIt = 50;
plan = [];

for n = 1:maxIt

  % take a random orientation
  ori = orientation.rand(cs);
  
  % simulate a corresponding noisy pattern
  pattern = det.simulatePattern(master,ori,1000,50);

  % approximate the pattern by a harmonic function -
  [pHarm,plan] = det.pattern2Fun(pattern,'bandwidth',128,'quadrature',plan,'keepPlan');

  % compute spherical convolution
  xcor = conv(masterHarm,pHarm,'bandwidth',bw);

  % correct for pattern shape
  xcor = FourierODF(xcor - sum(pHarm)*correction);

  % global search
  [~,id] = max(eval(xcor,S3G,'keepPlan'));
  %[~,id] = max(eval(xcor,S3G));

  % local search
  S3G_local = localGrid * S3G(id);
  [~,id] = max(eval(xcor,S3G_local));
    
  % store the misorientation angle
  omega = angle(S3G_local(id),ori)./degree;
  %disp([int2str(n) ' - ' xnum2str(omega)]);
  misorientation(n,1) = omega;
  
  progress(n,maxIt);
end

toc
eval(xcor.components{1},[],'killPlan')
S2FunHarmonic.quadrature([],'killPlan')

%save(['../sim/s2xcor_' int2str(bw) '_' int2str(100*S3G.resolution./degree) '_5.txt'],'misorientation','-ascii')


%% some plotting

close all
histogram(min(0.4,misorientation),15)
