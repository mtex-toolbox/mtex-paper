% MTEX script for the paper
%
%% Modelling and Testing Grain Boundary Plane Distributions
%
% by Ralf Hielscher, Rüdiger Kilian, Katharina Tinka Marquardt, Erik Wünsche
%
% * to run the script at least MTEX 6.1 has to be installed
% 
%% Section 3.1 - Kinetically driven boundary distributions

%% Initialize Neper

neper.init

neper.geometry = "cube(4,4,4)"

% neper.id = 529;

%% example 1 uniform ODF

cs = crystalSymmetry.load('quartz.cif');
odf = uniformODF(cs)

neper.morpho = "diameq:normal(10,2),aspratio(8,1,1)";

grains = neper.simulateGrains(4000,odf)

%save('grains31.mat','grains')

%%

load('grains31.mat')

%% visualize the 3d microstructure

close all

cKey = ipfColorKey(grains.CS);
cKey.inversePoleFigureDirection = zvector;

plot(grains,cKey.orientation2color(grains.meanOrientation),'micronbar','off')

pC = plottingConvention.default3D;
pC.setView
hx = text(2,0,0,"$x$",'FontSize',20,'Interpreter','latex');
hx.Position(3) = -0.1;
hy = text(0,2,0,"$y$",'FontSize',20,'Interpreter','latex');
hy.Position(3) = -0.15;

%exportgraphics(gca,'../pic/grainMacroUniform.png')
%exportgraphics(gca,'../pic/grainMacroTex.png')

%% the color Key

figure(2)
plot(cKey,'noTitle')

%exportgraphics(gca,'../pic/cKey.png')

%% the theoretical and simulated ODF

plotPDF(odf,Miller(0,0,1,cs))

setColorRange([0,8])

% exportgraphics(gca,'../pic/pfMacro.png')

%%

odfSim = calcDensity(grains.meanOrientation,'weights',grains.volume)

plotPDF(odfSim,Miller(0,0,1,cs))

mtexColorbar

setColorRange([0,8])

% exportgraphics(gca,'../pic/pfMacroSim.png')


%% the experimental specimen GBND

gB = grains.boundary('indexed')

gbndS = gB.calcGBND

%gbndS.how2plot = pC;

plot(gbndS)

%annotate(xvector,'label','X')
text([vector3d.X,vector3d.Y,vector3d.Z],{'X','Y','Z'},'BackgroundColor','w','tag','axesLabels')

mtexColorbar

%exportgraphics(gca,'../pic/gbndMacro.png')

%% the experimental crystal GBND

gbnd = gB.calcGBND(grains)

figure(1)
plot(gbnd,'complete','upper')
%mtexTitle('simulated GBND')
setColorRange([0.5,1.5])
mtexColorbar
%exportgraphics(gca,'../pic/bndAMacroUniformSim.png')

%% theoretical crystal GBND

gbndT = conv(odfSim,gbndS)

%nextAxis
plot(gbndT,'complete','upper')
%mtexTitle('theoretic GBND')
setColorRange([0.5,1.5])
%mtexColorbar
exportgraphics(gca,'../pic/bndAMacroUniform.png')


%% -------------------------------------------------------------
%% example 2 - textured microstructure

%define an ODF
odf = SO3FunHarmonic(SO3Fun.dubna);

odf = odf.symmetrise('SS',specimenSymmetry('orthorhombic'));

cs = odf.CS
%odf.SS.how2plot = pC;

plotPDF(odf,Miller({0,0,1},cs))
%exportgraphics(gca,'../pic/pdf001.png')

mtexColorbar

% update grain orientations
grains.meanOrientation = odf.discreteSample(length(grains))

grains = grains.update

%%

plot()


%% 

odfSim = calcDensity(grains.meanOrientation,'weights',grains.volume,'halfwidth',5*degree)

plotPDF(odfSim,Miller({0,0,1},cs))
%exportgraphics(gca,'../pic/pdf001.png')
setColorRange([0,8])
%mtexColorbar

% exportgraphics(gca,'../pic/pfMacroTexSim.png')

%%

plotPDF(odf,Miller({0,0,1},cs))

setColorRange([0,8])

%mtexColorbar

% exportgraphics(gca,'../pic/pfMacroTex.png')


%% simulated crystal GBND

gB = grains.boundary('indexed')

gbnd = gB.calcGBND(grains)

%figure(2)
plot(gbnd,'complete','upper')
%mtexTitle('simulated GBND')
setColorRange([0.5,1.5])
%mtexColorbar

%exportgraphics(gca,'../pic/bndAMacroTexSim.png')


%% theoretical crystal GBND

gbndT = conv(odf,gbndS)

%nextAxis
plot(gbndT,'complete','upper')
%mtexTitle('theoretic GBND')
setColorRange([0.5,1.5])
%mtexColorbar

%saveFigure('../pic/bndAMacroTex.png')


