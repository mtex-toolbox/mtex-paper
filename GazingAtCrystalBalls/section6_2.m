% Script file for the publication
%
%% Gazing at crystal balls  - 
%% electron backscatter diffraction indexing and cross correlation on the sphere
%
% R. Hielscher, F. Bartel, and T. B. Britton
%
% Ultramicroscopy, 2019
%
%% Section 6 Experimental Demonstration - The Cross Correlation Based Approach

%% Init Astro

Astro_FP='~/repo/AstroEBSD';
run([Astro_FP filesep 'start_AstroEBSD']);


%% Load and setup the the master pattern

BinFile='Fe-Iron-alpha.xml.BIN'; isHex=0;
cifname='Fe-Iron-alpha.cif';
[screen_int,facedata] = Cube_Generate(BinFile,isHex);

% as a function handle
fullPatternS2 = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,isHex));

% as hamonic expansion
cs = loadCIF(cifname);
patternFunHarm = S2FunHarmonicSym.quadrature(fullPatternS2,cs,'bandwidth',512);


%% load the h5 and set up background correction

InputUser.EBSD_File='data';
InputUser.HDF5_folder='';

[MapData,MicroscopeData,Phase_BrukerData,EBSD_DataInfo] = BCF_HDF5(InputUser);

%read the map & convert to area data
Data_InputMap = EBSD_Map(MapData,MicroscopeData);

% define correction setting for Kikuchi pattern
correctionSettings


%% parameters

% the harmonic bandwidth for the cross correlation function
bw = 48;

% the global search grid for peak detection
S3G = equispacedSO3Grid(cs,'resolution',1.5*degree);

% the local search grid for peak detection
localGrid = localOrientationGrid(specimenSymmetry,specimenSymmetry,...
  1.5*degree,'resolution',0.05*degree);

%% a class representing the detector parameters

det = detector(MicroscopeData.PatternWidth/4,MicroscopeData.PatternHeight/4,...
  MapData.DD(1),[MapData.PCX(1),MapData.PCY(1)]);

% define a cut off function for the detector
mask = det.S2CutOffMask(0.05);

% compute the auto correlation of the pattern with the cutoff function
correction = 1/sum(mask) * conv(patternFunHarm, mask,'bandwidth',bw) ;


%% whoole map analysis

% the detected orientations
xcf_h = NaN(Data_InputMap.xpts,Data_InputMap.ypts);
rot = orientation.nan(S3G.CS,size(xcf_h));
plan = [];

tic
for px = 1:Data_InputMap.xpts

  progress(px,Data_InputMap.xpts);

  for py = 1:Data_InputMap.ypts

    Pat_Number = Data_InputMap.PMap(py,px);
    pattern = bReadEBSP(EBSD_DataInfo,Pat_Number);
    pattern = EBSP_BGCor( pattern, Settings_Cor2 ).';

    % show pattern
    % imagesc(pattern)

    % approximate the pattern by a harmonic function -
    [pHarm,plan] = det.pattern2Fun(pattern,'bandwidth',128,'quadrature',plan);
    
    % compute spherical convolution
    xcor = conv(patternFunHarm,pHarm,'bandwidth',bw);

    % correct for pattern shape
    xcor = FourierODF(xcor - sum(pHarm)*correction);

    % the next command is the most demanding
    % speed depends mostly on the bandwidth bw
    % and the size of the search grid S3G
    %[xcf_h(px,py),id] = max(eval(xcor,S3G));
    
    % use this if the loop is not run in parallel - makes things 50% faster
    [xcf_h(px,py),id] = max(eval(xcor,S3G,'keepPlan')); 
    
    % local refinement
    S3G_local = localGrid * S3G(id);
    [~,id] = max(eval(xcor,S3G_local));
    rot(px,py) = S3G_local(id);
            
  end
end

toc

%save s2xcorMap_refined

% use this if the loop is not run in parallel
% this is to kill the NSOFT plan
eval(xcor,orientation.id(cs),'killPlan');

%% generate EBSD variable

% build the coordinate maps
prop.x = double(-Data_InputMap.XBeam_Map);
prop.y = double(Data_InputMap.YBeam_Map);

ebsd = EBSD(rot.', ones(size(rot)).',{'notIndexed',cs},'options',prop);
ebsd.prop.quality = xcf_h.';

% post rotate orientation data
ebsd = rotate(ebsd,rotation('axis',xvector,'angle',-MicroscopeData.TotalTilt),'keepXY');


%% plot ipf map

%set up cooordinate conventions
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');


colorKey = ipfHSVKey(cs);
colorKey.inversePoleFigureDirection = xvector;
color = colorKey.orientation2color(ebsd(cs.mineral).orientations);

figure(1); 
plot(ebsd(cs.mineral),color,'micronbar','off');
%hold on
%plot(ebsd,ebsd.prop.quality)
%mtexColorMap black2white
%hold off

%saveFigure('../pic/s2xcorrMap_refined.png')


%% plot mis2mean map

[grains,ebsd.grainId] = calcGrains(ebsd,'angle',5*degree);
colorKey = axisAngleColorKey;
colorKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;

figure(1); 
color = colorKey.orientation2color(ebsd(cs.mineral).orientations);
plot(ebsd(cs.mineral),color,'micronbar','off');

hold on
plot(grains.boundary,'linewidth',2)
hold off

%saveFigure('../pic/s2xcorrAxisAngle_refined.png')
