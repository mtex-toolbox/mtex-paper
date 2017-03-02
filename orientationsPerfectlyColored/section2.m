% Script file for Section 2 of the publication
%
%% IPF coloring of crystal orientation data
%
% Gert Nolze and Ralf Hielscher 
%

%% Section 2.2 No Symmetry P = 1

%%
% Figure 3:
%
% Colorization of the unit sphere using fully saturated colors, i.e., S = 1
% with hue H being the azimuth of a point on the sphere and the lightning L
% being the polar angle.

cs = crystalSymmetry('1')
oM = ipdfHSVOrientationMapping(cs); 
%oM.colorPostRotation = rotation('Euler',-pi/2,pi/2,0,'ZYZ');
oM.grayValue = 0.2;
oM.grayGradient = 0.2;
plot(oM,'noTitle','grid','grid_res',30*degree)

mtexFig = gcm;
sP = getappdata(mtexFig.children(1),'sphericalPlot');
set(sP.grid,'color',[0 0 0]);
sP = getappdata(mtexFig.children(2),'sphericalPlot');
set(sP.grid,'color',[0 0 0]);

%saveFigure('../pic/1.png')



%% Section 2.4 Colorizing spherical triangles and lunes
%
% Figure 4a

cs = crystalSymmetry('mmm');
sR = cs.fundamentalSector;

oM = ipdfHSVOrientationMapping(cs)
plot(oM,'noTitle','autoAlignText','Marker','none')

plotzOutOfPlane
v = plotS2Grid(sR,'resolution',0.5*degree);

[r,rho] = sR.polarCoordinates(v,sR.center);

hold on
contour(v,r,'contours',7)
hold on
contour(v,rho,'contours',10)
hold on
plot(sR,'linewidth',2)


plot(sR.center,'label','p','backgroundcolor',...
  'w','marker','s','markercolor','k')
hold off
%saveFigure('../pic/polarTrianglemmm.png')

%%
% Figure 4b

cs = crystalSymmetry('m-3m');
sR = cs.fundamentalSector;

oM = ipdfHSVOrientationMapping(cs)

plot(oM,'noTitle','autoAlignText','Marker','none')

v = plotS2Grid(sR,'resolution',0.5*degree);

[r,rho] = sR.polarCoordinates(v,sR.center);

hold on
contour(v,r,'contours',7)
hold on
contour(v,rho,'contours',10)
hold on
plot(sR,'linewidth',2)


plot(sR.center,'label','p','backgroundcolor','w','marker','s','markercolor','k')
hold off
%saveFigure('../pic/polarTrianglem-3m.png')

%%
% Figure 4c

cs = crystalSymmetry('31m');
sR = cs.fundamentalSector;

oM = ipdfHSVOrientationMapping(cs)

plot(oM,'noTitle','autoAlignText','Marker','none')

v = plotS2Grid(sR,'resolution',0.5*degree);

[r,rho] = sR.polarCoordinates(v,sR.center);

hold on
contour(v,r,'contours',7)
hold on
contour(v,rho,'contours',13)
hold on
plot(sR,'linewidth',2)

plot(sR.center,'label','p','backgroundcolor','w','marker','s','markercolor','k')
hold off
%saveFigure('../pic/polarTriangle3m.png')

%% Section 2.6 Spherical triangles without reflection boundaries
%
% Figure 5a

% define a ipf key for m-3
cs = crystalSymmetry('m-3')
oM = ipdfHSVOrientationMapping(cs);

% without additional reflections
oM.refl = [];

% and with respect to the true fundamental sector
oM.sR = cs.fundamentalSector;

% and with respect to the true center
oM.whiteCenter = oM.sR.center;

% plot the ipf key
plot(oM,'complete','upper','noLabel','resolution',0.25*degree)
hold on
plot(cs)
hold on
plot(cs.fundamentalSector,'color','w','linewidth',3)
hold off

%saveFigure('../pic/m-3A.png')

%%
%
% Figure 5b

% use a higher symmetric key
cs2 = crystalSymmetry('m-3m')
oM = ipdfHSVOrientationMapping(cs2);
plot(oM,'complete','upper','noLabel','resolution',1*degree)
hold on
plot(cs)
hold on
plot(cs.fundamentalSector,'color','w','linewidth',3)
hold off
%saveFigure('../pic/m-3B.png')

%%
%
% Figure 5c

% the correct key with additional reflections
oM = ipdfHSVOrientationMapping(cs);
plot(oM,'complete','upper','noLabel','resolution',1*degree)
hold on
plot(cs)
hold on
plot(cs.fundamentalSector,'color','w','linewidth',2)
hold off
%saveFigure('../pic/m-3C.png')

%% Section 2.7 The exceptions: -1, -3, -4
%
%%
% Figure 6a - ipf key for -1 with color discontinuities

cs = crystalSymmetry('-1')
oM = ipdfHSVOrientationMapping(cs)

plot(oM,'complete','noLabel','notitle')

% saveFigure('../full/-1a.png')

%%

% Figure 6b - ambiguis ipf key for -1 without color discontinuities

cs = crystalSymmetry('-1')
oM = ipdfHSVOrientationMapping(cs)
oM.maxAngle = pi/2
oM.maxAngle = pi/2
plot(oM,'complete','noLabel','notitle')

% saveFigure('../full/-1b.png')

%%

% Figure 7a - ipf key for -3 with color discontinuities

cs = crystalSymmetry('-3','a||x')
oM = ipdfHSVOrientationMapping(cs)

plot(oM,'complete','nolabel','notitle')
hold on
plot(cs)
hold off

% saveFigure('../full/-3a.png')

%%

% Figure 7b - ipf key for -4 with color discontinuities

cs = crystalSymmetry('-4')
oM = ipdfHSVOrientationMapping(cs)

plot(oM,'complete','noLabel','notitle')
hold on
plot(cs)
hold off

% saveFigure('../full/-4a.png')
