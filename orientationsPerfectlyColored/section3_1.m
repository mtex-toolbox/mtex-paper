% Script file for Section 3.1 of the publication
%
%% IPF coloring of crystal orientation data
%
% Gert Nolze and Ralf Hielscher 
%

%% Data import
% import ideal Kurdjumov-Sachs misorientations

plotx2east
ebsd = EBSD.load('K-SX.ctf','convertspatial2EulerReferenceFrame')
ori = ebsd.orientations;

grains = calcGrains(ebsd,'unitCell');



%% Figure 8
%
% Figure 8a
% plot the misorientations in an EBSD map and colorize them according to
% different reference directions - first with respect to X


ipfKey = ipfHSVKey(ebsd.CS);

% here the reference direction is set
ipfKey.inversePoleFigureDirection = vector3d(1,0,0);
ipfKey.grayGradient = 0.25;
ipfKey.grayValue = 0.2;

plot(ebsd,ipfKey.orientation2color(ori),'EdgeColor','k','micronbar',false,'linewidth',2)

id = [2,8,5,11];
ori = grains.meanOrientation;
ind = ismember(1:24,id);

hold on
%plot(grains(id).boundary,'linecolor','b','linewidth',5)

% annotate grain ids
for i = 1:4
  [x,y] = grains(id(i)).centroid;
  text(x,y,xnum2str(i),'HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',50)
end
hold off

%saveFigure('../pic/K-SX-IPFX.png')

%%
% Figure 8b
% next colorize the orientations with respect to Y

ipfKey.inversePoleFigureDirection = yvector;
plot(ebsd,ipfKey.orientation2color(ori),'EdgeColor','k','micronbar',false,'linewidth',2)

hold on
%plot(grains(id).boundary,'linecolor','b','linewidth',5)

for i = 1:4
  [x,y] = grains(id(i)).centroid;
  text(x,y,xnum2str(i),'HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',50)
end
hold off

%saveFigure('../pic/K-SX-IPFY.png')

%%
% Figure 8c
% and finally with respect to z

ipfKey.inversePoleFigureDirection = zvector;
plot(ebsd,ipfKey.orientation2color(ori),'EdgeColor','k','micronbar',false,'linewidth',2)

hold on
%plot(grains(id).boundary,'linecolor','b','linewidth',5)

for i = 1:4
  [x,y] = grains(id(i)).centroid;
  text(x,y,xnum2str(i),'HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',50)
end
hold off

%saveFigure('../pic/K-SX-IPFZ.png')


%%
% Figure 8d
%
% plot the colorized orientations into a pole figure

ipfKey.inversePoleFigureDirection = yvector;

plotPDF(ori,ipfKey.orientation2color(ori),Miller(1,0,0,ebsd.CS),...
  'projection','eangle','noTitle','MarkerEdgeColor',0.6*[1 1 1],...
  'MarkerSize',7,'grid','on')

%hold on
%circle([xvector,yvector],'linecolor',0.3*[1 1 1])

text([xvector,yvector],{'X','Y'},'backgroundcolor','w')

% the following lines annotates the orientations
%hold on
%plotPDF(ori(ind),Miller(1,0,0,ebsd.CS),...
%  'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5,'noTitle')
%plotPDF(ori(ind),oM.orientation2color(ori(ind)),Miller(1,0,0,ebsd.CS),...
%  'MarkerEdgeColor','b','MarkerSize',4,'noTitle')
%
%hold on

hold off

drawNow(gcm,'figSize','large')

% saveFigure('../pic/K-SX-PF100.pdf')

%%
% Figure 8d
% mark the orientations in the IPF-Z key


ori.CS = ipfKey.CS1;
ipfKey.inversePoleFigureDirection = zvector;
plot(ipfKey,'noTitle','autoAlignText','Marker','none')
hold on
plot(ori,'markerFaceColor','k','MarkerEdgeColor','k','marker','o','MarkerSize',7)
hold off

%saveFigure('../pic/K-SX-IPFKEY.png')


%% Figure 9
% Figure 9a, 9b, 9c

% rotate the orientations a bit
ori_rot = rotation('Euler',25*degree,15*degree,55*degree) * ori;

% change to xvector, yvector, zvector to generate a, b, c
ipfKey.inversePoleFigureDirection = xvector;

plot(ebsd,ipfKey.orientation2color(ori_rot),'EdgeColor','k','micronbar',false,'linewidth',2)

% saveFigure('../pic/K-SX-ROT-IPFX.png')

%%
% Figure 9d

plotPDF(ori_rot,ipfKey.orientation2color(ori_rot),Miller(1,0,0,ebsd.CS),...
  'projection','eangle','noTitle','grid','on',...
  'MarkerEdgeColor',0.6*[1 1 1],'MarkerSize',7)

text([xvector,yvector],{'X','Y'},'backgroundcolor','w')

drawNow(gcm,'figSize','large')

% saveFigure('../pic/K-SX-ROT-PF100.pdf')


%%
% Figure 9e

ipfKey.inversePoleFigureDirection = yvector;
plot(ipfKey,'noTitle','autoAlignText','Marker','none')
hold on
plot(ori_rot,'markerFaceColor','k','MarkerEdgeColor','k','marker','o','MarkerSize',7)
hold off

%saveFigure('../pic/K-SX-IPFROTKEY.png')
