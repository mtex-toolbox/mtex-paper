% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%

%% Section 4.4 Orientation Determination

%% init Astro

Astro_FP='~/repo/AstroEBSD';
run([Astro_FP filesep 'start_AstroEBSD']);

[ UCell,Crystal_Family,Crystal_LUT,Settings_LUT] = ...
  Phase_Builder({'Ferrite'},[Astro_FP filesep 'phases'] );

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

% define a cut off function for the detector
mask = det.S2CutOffMask(0);
maskHarm = S2FunHarmonic.quadrature(mask);


%% define a band profile

% profile 1
profile = @(x) 1.5*exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-86*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-94*degree).^2./(2*degree).^2);

profileHarm = S2Kernel.quadrature(profile);


%% Loop for multiple runs

% define the harmonic degree
bw = 64;
setMTEXpref('maxBandwidth',bw);
numBands = 12;

% 
correction = max(maskHarm.radon,0.1);

tic
%profile on
%maxIt = 5000;
maxIt = 500;
plan = [];

for n = 1:maxIt

  % simulate a pattern
  ori = orientation.rand(cs);
  pattern = det.simulatePattern(master,ori,1000,50);

  % approximate the pattern by a harmonic function -
  [pHarm,plan] = det.pattern2Fun(pattern,'bandwidth',bw*2,'quadrature',plan);

  % corrected Radon transform
  %RadonPHarm= pHarm.radon ./ correction;
  RadonPHarm= conv(pHarm,profileHarm) ./ correction;

  % Peak detection
  peaks = S2PeakDetection(RadonPHarm,det,numBands);

  % Use Astro for indexing
  [rot,bands] = EBSP_Index(squeeze(double(peaks)), ...
    Crystal_LUT{1}, Settings_LUT{1}.thresh_trig, UCell{1});

  ori_found = orientation('Euler',rot.eang,cs);
  omega = angle(ori,ori_found) ./ degree;
  %if omega > 10, remember = ori; end %#ok<BDSCI>
  %if omega > 3, break; end %#ok<BDSCI>
  misorientation(n,1) = omega; %#ok<SAGROW>
    
  progress(n,maxIt);
end

%profview
toc

%save(['../sim/band_' int2str(bw) '_' int2str(numBands) '.txt'],'misorientation','-ascii')

%% misorientation angle diagramm -> Figure 6
clf
histogram(min(misorientation(1:n-1),1),30)
xlabel('misorientation angle in degree')
set(gca,'FontSize',16)
% saveFigure('../pic/bandHist128.pdf')

%save bandDetectionSim