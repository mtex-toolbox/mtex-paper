% Script file for Section 1 of the publication
%
%% IPF coloring of crystal orientation data
%
% Gert Nolze and Ralf Hielscher 
%

%% Section 1.1
%
% Import the quartz data and consider only the largest grain

% load the ebsd data
ebsd = loadEBSD('Quartz-175B_05.ctf','convertEuler2SpatialReferenceFrame');

% restrict to the quartz phase
ebsd = ebsd('Quartz');

% compute grains
[grains,ebsd.grainId] = calcGrains(ebsd)

% and select the big one
[~,id] = max(grains.grainSize);
ebsd = ebsd(grains(id));
ori = ebsd.orientations;

% some plotting conventions
plotzOutOfPlane
plotx2east


%% Figure 1a
%
% Euler RGB color map of the quartz grain

% define the orientation to color map
oM = BungeRGBOrientationMapping(ori);

% rotate the data to make the color jump better visible
ebsd_rot = rotate(ebsd,rotation('Euler',90*degree,180*degree,0),'keepXY');

% plot the orientation map colorized by Bunge Euler angles
plot(ebsd,oM.orientation2color(ebsd_rot.orientations))
%saveFigure('../quartz/quartz_RGB.png')


%% Figure 1b
%
% indiviudal orientations of the quartz grain as phi1 - phi2 plot

close all
% reduce number of orientations to be plotted to 1000 randomly choosen one
ori2 = discreteSample(ebsd_rot.orientations,1000);

% project Euler angles to the fundamental sector
% phi2 = [0..120], Phi = [0..90] 
[phi1,Phi,phi2] = project2EulerFR(ori2,'Bunge');

% scatter plot phi1 - phi2 colorized by the RGB as the map from Fig. 1a
scatter(phi1./degree, phi2./degree,10, oM.orientation2color(ori2))

% make the plot nnice
xlabel('$\varphi_1$','fontSize',20,'interpreter','LaTeX')
zlabel('\Phi')
ylabel('$\varphi_2$','fontSize',20,'interpreter','LaTeX')
set(gca,'FontSize',20,'xtick',[0 120 240 360],'ytick',[0 40 80 120])
xlim([0,360])
ylim([0,120])
box on

%saveFigure('../quartz/euler.pdf')

%% Section 1.2

%% Figure 1c
% Colorize the quartz grain according to the TSL ipf key

% define the TSL IPF key
oM = TSLOrientationMapping(ebsd);

% set the reference direction
oM.inversePoleFigureDirection = yvector;

% plot the orientation map
plotzIntoPlane
plot(ebsd,oM.orientation2color(ori))
%saveFigure('../quartz/quartz_TSL.png')

%% Figur 1d
% Plot the TSL IPF key

% plot the key
plotzOutOfPlane
plot(oM,'noTitle','autoAlignText','Marker','none')

% plot the orientations on top of the key
hold on
plot(ori,'MarkerFaceColor','none','MarkerEdgeColor','r','points',20)
hold off

%saveFigure('../quartz/quartz_OM_TSL.png')

%% Section 1.3

csNames = {'-1','112/m','mmm','4/m','4/mmm','-3','-3m','6/m','6/mmm','m-3'}

cs = crystalSymmetry(csNames{end-1})
oM = TSLOrientationMapping(cs)

plot(oM,'3d')
%plot(oM,'upper','complete','noLabel','noTitle')
%hold on
%plot(cs)
%hold off
