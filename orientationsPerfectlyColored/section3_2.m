% Script file for Section 3.2 of the publication
%
%% IPF coloring of crystal orientation data
%
% Gert Nolze and Ralf Hielscher 
%


%% Section 3.2.1 - Encoding for centrosymmetric point groups.

% data import
plotx2east
plotzOutOfPlane
ebsd = loadEBSD('K-SX.ctf','convertspatial2EulerReferenceFrame')
ori = ebsd.orientations;

grains = calcGrains(ebsd,'unitCell');

%%
% Figure 10a

oM = ipdfHSVOrientationMapping(ebsd.CS.properGroup);
oM.inversePoleFigureDirection = xvector;
oM.grayGradient = 0.25;
oM.grayValue = 0.2;

plot(ebsd,oM.orientation2color(ori),'EdgeColor','k','micronbar',false,'linewidth',2)

id = [2,8,5,11];
ori = grains.meanOrientation;
ind = ismember(1:24,id);

hold on
for i = 1:4
  [x,y] = grains(id(i)).centroid;
  text(x,y,xnum2str(i),'HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',50)
end
hold off

%saveFigure('../pic/K-SX-IPFX_proper.png')

%%
% Figure 10b

oM.inversePoleFigureDirection = yvector;
plot(ebsd,oM.orientation2color(ori),'EdgeColor','k','micronbar',false,'linewidth',2)

hold on
for i = 1:4
  [x,y] = grains(id(i)).centroid;
  text(x,y,xnum2str(i),'HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',50)
end
hold off

%saveFigure('../pic/K-SX-IPFY_proper.png')

%%
% Figure 10c

oM.inversePoleFigureDirection = zvector;
plot(ebsd,oM.orientation2color(ori),'EdgeColor','k','micronbar',false,'linewidth',2)

hold on
for i = 1:4
  [x,y] = grains(id(i)).centroid;
  text(x,y,xnum2str(i),'HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',50)
end
hold off

%saveFigure('../pic/K-SX-IPFZ_proper.png')

%%
% Figure 10d

oM.inversePoleFigureDirection = zvector;

plotPDF(ori,oM.orientation2color(ori),Miller(1,0,0,ebsd.CS),...
  'projection','eangle','noTitle','MarkerEdgeColor',0.6*[1 1 1],...
  'MarkerSize',7,'grid','on')

%hold on
%circle([xvector,yvector],'linecolor',0.3*[1 1 1])

text([xvector,yvector],{'X','Y'},'backgroundcolor','w')

%hold on
%plotPDF(ori(ind),Miller(1,0,0,ebsd.CS),...
%  'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5,'noTitle')
%plotPDF(ori(ind),oM.orientation2color(ori(ind)),Miller(1,0,0,ebsd.CS),...
%  'MarkerEdgeColor','b','MarkerSize',4,'noTitle')
%
%hold on

hold off

drawNow(gcm,'figSize','large')

%saveFigure('../pic/K-SX-PF100_proper.pdf')

%%
% Figure 10e

ori.CS = oM.CS1;
oM.inversePoleFigureDirection = zvector
plot(oM,'noTitle','autoAlignText','Marker','none')
hold on
plot(ori,'markerFaceColor','k','MarkerEdgeColor','k','marker','o','MarkerSize',7)
hold off

%saveFigure('../pic/K-SX-IPFKEY_proper.png')


%% Section 3.2.2 Encoding for enantiomorphic point groups

%%
% data peperation

% load the ebsd data
ebsd = loadEBSD('Quartz-175B_05.ctf','convertEuler2SpatialReferenceFrame');

% restrict to the quartz phase
ebsd = ebsd('Quartz');

% compute grains
[grains,ebsd.grainId] = calcGrains(ebsd);

% select only the big grain and the inclusion
grains = grains([18,74]);
ebsd = ebsd(grains);
grains = calcGrains(ebsd);


%%
% Figure 11a

oM = ipdfHSVOrientationMapping(ebsd.CS);
oM.grayValue = [0.2,0.5];
oM.grayGradient = 0.5

plotzOutOfPlane
plot(oM,'noTitle','autoAlignText','Marker','none')

ori = grains.meanOrientation;

plotzOutOfPlane
plot(oM,'noTitle','autoAlignText','Marker','none')

args = {'MarkerSize',10,'symmetrised'};

hold on
plot(ori(1) .\ vector3d(0,0,1),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k','marker','d')
plot(ori(2) .\ vector3d(0,0,1),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w','marker','d')

plot(ori(1) .\ vector3d(1,0,0),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k','marker','s')
plot(ori(2) .\ vector3d(1,0,0),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w','marker','s')

plot(ori(1) .\ vector3d(0,1,1),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(ori(2) .\ vector3d(0,1,1),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w')

plot(ori(1) .\ vector3d(1,1,1),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k','marker','^')
plot(ori(2) .\ vector3d(1,1,1),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w','marker','^')
hold off

%saveFigure('../quartz/quartz_OM_Laue.png')

%%

% Figure 11c
plot(oM,'noTitle','complete','upper')
hold on
plot(oM.CS1)
hold off

%saveFigure('../quartz/-3m1.png')

%%
% Figure 11b

oM = ipdfHSVOrientationMapping(ebsd.CS.properGroup);
oM.grayValue = [0.2,0.5];
oM.grayGradient = 0.25;

ori = grains.meanOrientation;
ori.CS = ori.CS.properGroup;

plotzOutOfPlane
plot(oM,'noTitle','autoAlignText','Marker','none')

args = {'MarkerSize',10,'symmetrised'};

hold on
plot(ori(1) .\ vector3d(0,0,1),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k','marker','d')
plot(ori(2) .\ vector3d(0,0,1),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w','marker','d')

plot(ori(1) .\ vector3d(1,0,0),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k','marker','s')
plot(ori(2) .\ vector3d(1,0,0),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w','marker','s')

plot(ori(1) .\ vector3d(0,1,1),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(ori(2) .\ vector3d(0,1,1),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w')

plot(ori(1) .\ vector3d(1,1,1),args{:},'MarkerFaceColor','w','MarkerEdgeColor','k','marker','^')
plot(ori(2) .\ vector3d(1,1,1),args{:},'MarkerFaceColor','k','MarkerEdgeColor','w','marker','^')
hold off

%saveFigure('../quartz/quartz_OM_proper.png')

%%

% Figure 11d
plot(oM,'noTitle','complete','upper')
hold on
plot(oM.CS1)
hold off

%saveFigure('../quartz/321.png')

%%
% Figure 12a

r = xvector;
% r = zvector;
% r = yvector + zvector;
% r = xvector + yvector + zvector;

oM = ipdfHSVOrientationMapping(ori.CS.Laue);
oM.inversePoleFigureDirection = r;
oM.grayValue = [0.2,0.5];
oM.grayGradient = 0.5;
plotzIntoPlane
plot(ebsd,oM.orientation2color(ebsd.orientations))
hold on
plot(grains(1).boundary,'linewidth',1.5)
hold off

%saveFigure('../quartz/quartzTwin100.png')
%saveFigure('../quartz/quartzTwin001.png')
%saveFigure('../quartz/quartzTwin011.png')
%saveFigure('../quartz/quartzTwin111.png')

%%
% Figure 12b

r = xvector;
% r = zvector;
% r = yvector + zvector;
% r = xvector + yvector + zvector;

oM = ipdfHSVOrientationMapping(ori.CS.properGroup);
oM.inversePoleFigureDirection = r;
plotzIntoPlane
plot(ebsd,oM.orientation2color(ebsd.orientations))
hold on
plot(grains(1).boundary,'linewidth',1.5)
hold off

%saveFigure('../quartz/quartzTwin100Proper.png')
%saveFigure('../quartz/quartzTwin001Proper.png')
%saveFigure('../quartz/quartzTwin011Proper.png')
%saveFigure('../quartz/quartzTwin111Proper.png')

%% Section 3.2.3 - Encoding for neither centrosymmetric  nor enantiomorphic point groups

% Import the Data

% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

CS = crystalSymmetry('-43m', [5.4505 5.4505 5.4505], 'mineral', 'GaP');
ebsd = loadEBSD_generic('GaP320_XCDS.BEST.MTXT','CS',CS,...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'phase' 'x' 'y' 'xc' 'int' 'kik'})

grains = calcGrains(ebsd,'angle',5*degree);

% identify three grains
ids = [6218,21411,27615];
grainIds = grains.findByLocation([ebsd(ids).x,ebsd(ids).y]);

%%
% Figure 13a

% consider rotational part of the Laue group of 
oM = ipdfHSVOrientationMapping(CS.properGroup);
oM.inversePoleFigureDirection = vector3d.X;
oM.grayGradient = 0.25;
oM.grayValue = [0.2 0.5];

figure(1)
plot(ebsd,oM.orientation2color(ebsd.orientations),'micronBar',false)
hold on
plot(grains.boundary)
hold off

% mark some grains
c = grains(grainIds).centroid;
text(c(:,1),c(:,2),{'A','B','C'},'fontSize',10,'backgroundcolor','w')

%saveFigure('../GaP/432map.png')

%%
% Figure 14a
% this is the colored fundamental sector
figure(2)
cs = crystalSymmetry('23', [5.4505 5.4505 5.4505], 'mineral', 'GaP');
sR = cs.fundamentalSector;
plot(oM,sR)
hold on
ori = orientation(grains(grainIds).meanOrientation,CS.properSubGroup);
plot(ori([1,3]) \ xvector,'symmetrised','label',{'A','C'},...
  'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','w','backgroundColor','w')
plot(ori(2) \ xvector,'symmetrised','label',{'B'},'textaboveMarker',...
  'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','w','backgroundColor','w')
hold off

%saveFigure('../GaP/432key.png')

%%
% Figure 13b

% consider rotational part of the point group
oM = ipdfHSVOrientationMapping(CS.properSubGroup);
oM.inversePoleFigureDirection = vector3d.X;
oM.grayGradient = 0.25;
oM.grayValue = [0.2 0.5];

figure(3)
plot(ebsd,oM.orientation2color(ebsd.orientations),'micronBar',false)
hold on
plot(grains.boundary)
hold off

% mark some grains
c = grains(grainIds).centroid;
text(c(:,1),c(:,2),{'A','B','C'},'backgroundcolor','w')

%saveFigure('../GaP/23map.png')

%%
% Figure 14b

% this is the colored fundamental sector
figure(4)
plot(oM)
hold on
ori = orientation(grains(grainIds).meanOrientation,CS.properSubGroup);
plot(ori([1,3]) \ xvector,'symmetrised','label',{'A','C'},...
  'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','w','backgroundColor','w')
plot(ori(2) \ xvector,'symmetrised','label',{'B'},'textaboveMarker',...
  'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','w','backgroundColor','w')
hold off

%saveFigure('../GaP/23key.png')

%%

angle(grains(grainIds(1)).meanOrientation,grains(grainIds(2)).meanOrientation) ./ degree
angle(grains(grainIds(2)).meanOrientation,grains(grainIds(3)).meanOrientation) ./ degree
angle(grains(grainIds(1)).meanOrientation,grains(grainIds(3)).meanOrientation) ./ degree

