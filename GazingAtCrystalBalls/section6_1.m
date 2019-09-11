% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%
%% Section 6 Experimental Demonstration - The Radon Transform Based Approach

clear; close all; home;
tic

%set up cooordinate conventions
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

%% init Astro

cs = loadCIF('Fe-Iron-alpha.cif');

Astro_FP='~/repo/AstroEBSD';
run([Astro_FP filesep 'start_AstroEBSD']);

[ UCell,Crystal_Family,Crystal_LUT,Settings_LUT] = ...
  Phase_Builder({'Ferrite'},[Astro_FP filesep 'phases'] );

%% set up the detector

InputUser.EBSD_File='data';
InputUser.HDF5_folder='';

[MapData,MicroscopeData,Phase_BrukerData,EBSD_DataInfo] = BCF_HDF5( InputUser );

det = detector(MicroscopeData.PatternWidth/4,MicroscopeData.PatternHeight/4,...
  MapData.DD(1),[MapData.PCX(1),MapData.PCY(1)]);

% define a cut off function for the detector
mask = det.S2CutOffMask(0);
maskHarm = S2FunHarmonic.quadrature(mask);

% define correction setting for Kikuchi pattern
correctionSettings

%% read the map & convert to area data

Data_InputMap = EBSD_Map(MapData,MicroscopeData);

%% define a band profile for band detection

profile = @(x) exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-87*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-93*degree).^2./(2*degree).^2);

%profile = @(x) exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
%  exp(-(acos(x)-86*degree).^2./(2*degree).^2) - ...
%  exp(-(acos(x)-94*degree).^2./(2*degree).^2);

%f = @(x) exp(-(acos(x)-90*degree).^2./(5*degree).^2) ;

% express the profile in Legendre coefficients
profileHarm = S2Kernel.quadrature(profile);

%% whoole map analysis
setMTEXpref('maxBandwidth',64)
numBands = 10;
correction = max(maskHarm.radon,0.1);

ori = orientation.nan(Data_InputMap.xpts,Data_InputMap.ypts,cs);
plan = [];

for px = 1:Data_InputMap.xpts

  progress(px,Data_InputMap.xpts);

  for py = 1:Data_InputMap.ypts

    Pat_Number = Data_InputMap.PMap(py,px);
    pattern = bReadEBSP(EBSD_DataInfo,Pat_Number);
    pattern = EBSP_BGCor( pattern, Settings_Cor2 ).';
    
    % approximate the pattern by a harmonic function -
    [pHarm,plan] = det.pattern2Fun(pattern,'bandwidth',128,'quadrature',plan);

    % corrected Radon transform
    RadonPHarm= conv(pHarm,profileHarm) ./ correction;

    % peak detection - last argument is the number of peaks to detect
    peaks = S2PeakDetection(RadonPHarm,det,numBands);
    
    % Use Astro for indexing
    [rot,bands] = EBSP_Index(squeeze(double(peaks)), ...
      Crystal_LUT{1}, Settings_LUT{1}.thresh_trig, UCell{1});

    ori(px,py) = orientation('Euler',rot.eang,cs);
            
  end
end

% save bandMap64_10

%% generate EBSD variable

% build the coordinate maps
prop.x = double(-Data_InputMap.XBeam_Map); %have to transpose - not sure why...
prop.y = double(Data_InputMap.YBeam_Map);

ebsd = EBSD(ori.', ones(size(ori)),{'notIndexed',cs},'options',prop);

% post rotate orientation data
ebsd = rotate(ebsd,rotation('axis',xvector,'angle',-MicroscopeData.TotalTilt),'keepXY');


%% plotting as ipf map

%colorKey = ipfTSLKey(cs);
colorKey = ipfHSVKey(cs);
colorKey.inversePoleFigureDirection = xvector;
color = colorKey.orientation2color(ebsd('indexed').orientations);

%figure(1); 
plot(ebsd(cs.mineral),color,'micronbar','off');
%saveFigure('../pic/mapBandDetection.png')

%% plotting as axis angle map with respect to grain mean orientation

[grains,ebsd.grainId] = calcGrains(ebsd,'angle',5*degree);
colorKey = axisAngleColorKey;
colorKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;

figure(1); 
color = colorKey.orientation2color(ebsd(cs.mineral).orientations);
plot(ebsd(cs.mineral),color,'micronbar','off');

hold on
plot(grains.boundary,'linewidth',2)
hold off

%saveFigure('bandDetectionAxisAngle.png')

%% Bruker indexed data for comparison as ipf map

figure(2)

prop.x = double(-Data_InputMap.XBeam_Map); %have to transpose - not sure why...
prop.y = double(Data_InputMap.YBeam_Map);

ebsd2 = EBSD(rotation('Euler',Data_InputMap.phi1*degree,Data_InputMap.PHI*degree,Data_InputMap.phi2*degree),...
    ones(size(Data_InputMap.phi1)),{'notIndexed',cs},'options',prop);

% ebsd2 = ebsd;
% ebsd2.rotations = ...
%   rotation('Euler',Data_InputMap.phi1*degree,Data_InputMap.PHI*degree,Data_InputMap.phi2*degree);

colorKey = ipfHSVKey(cs);
colorKey.inversePoleFigureDirection = xvector;

plot(ebsd2,colorKey.orientation2color(ebsd2.orientations),'micronbar','off')

%% as axis angle map with respect to grain mean orientation

colorKey = axisAngleColorKey;
[grains2,ebsd2.grainId] = calcGrains(ebsd2);
colorKey.oriRef = grains2(ebsd2('indexed').grainId).meanOrientation;


figure(2); 
color = colorKey.orientation2color(ebsd2(cs.mineral).orientations);
plot(ebsd2(cs.mineral),color,'micronbar','off');

hold on
plot(grains2.boundary)
hold off

%saveFigure('BrukerAxisAngle.png')


%% peak visualization

figure(4)
plot(RadonPHarm,'pcolor','resolution',0.25*degree,'upper')
%colormap gray
annotate(peaks,'markerSize',10,'antipodal','MarkerEdgeColor','k','MarkerFaceColor','none')


%% band visualization

figure(1)
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap gray

circle(peaks,'color','r')
