% MTEX script for the paper
%
%% Modeling and Testing Grain Boundary Plane Distributions
%
% by Ralf Hielscher, Rüdiger Kilian, Katharina Tinka Marquardt, Erik Wünsche
%
% * to run the script at least MTEX 6.1 has to be installed
% 
%% Section 3.2 - Crystallographically Driven Boundary Distributions
% 

%% Example 1 - no texture

cs = crystalSymmetry('432');

% Sigma 3 twinning
moriRef = orientation.byAxisAngle(Miller(1,1,1,cs),180*degree);

% the habit theoretic plane distribution
bndA = 0.4 * S2FunRBF(Miller(1,1,1,cs),'halfwidth',2.5*degree) ...
  + 0.6*S2FunRBF(Miller(1,0,1,cs),'halfwidth',5*degree);
bndA.antipodal = true;
bndA.CS = cs;

plot(bndA)
mtexColorbar
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)

%exportgraphics(gca,'../pic/bndADynamic.png')

%% the effective GBCD

% compute intertwining symmetry elements
[~,~,csRed] = project2FundamentalRegion(moriRef,moriRef,moriRef);

bndB = rotate(bndA,inv(moriRef));
bndAB = 0.5*(bndA + bndB);

plot(bndA.symmetrise(csRed),'upper')
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)
% exportgraphics(gca,'../pic/bndADynamic.png')

%%

plot(bndB.symmetrise(csRed),'upper')
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)
mtexColorbar

% exportgraphics(gca,'../pic/bndBDynamic.png')

%%

plot(bndAB.symmetrise(csRed),'upper')

%%

plot(bndAB.symmetrise(csRed),'upper')
% exportgraphics(gca,'../pic/bndABDynamic.png')


%% simulate microstructure

neper.morpho = 'diameq:normal(10,1),aspratio(1,1,1)';
grains = neper.simulateTwinGrains(5000,moriRef,uniformODF(cs),bndA)

% save("grains32a.mat",'grains','moriRef','bndA','cs')

%%

load("grains32a.mat")


%%

plot(grains,grains.meanOrientation)

pC = plottingConvention.default3D;
pC.setView

%exportgraphics(gca,'../pic/grainDynamicUniform.png')

%% experimental GBCD
% since mori = inv(mori) we can not distinguish between gbndA and gbndB

gB = grains.boundary('indexed');

gbndAB = gB.calcGBND(grains,moriRef,'halfwidth',1.5*degree)

figure(1)
plot(gbndAB,'complete','upper')
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)
mtexColorbar
%mtexTitle('simulated crystalographic GBND')

% exportgraphics(gca,'../pic/gbcdDynamic.png')

%% ODF estimation

odfA = uniformODF(cs)
odfEst = calcDensity(grains.meanOrientation,'weights',grains.volume)

%%

plotPDF(odfEst,Miller(1,0,0,cs))
setColorRange([0,3.3])
mtexColorbar

% exportgraphics(gca,'../pic/pfDynamicSim.png')

%%

plotPDF(uniformODF(cs),Miller(1,0,0,cs))
setColorRange([0,3.3])

% exportgraphics(gca,'../pic/pfDynamic.png')

%% experimental specimen GBND

gB = grains.boundary('indexed');
gbndS = gB.calcGBND('halfwidth',10*degree)

plot(gbndS,'figSize','normal')
mtexColorbar
setColorRange([0,4.2])
%mtexTitle('experimental specimen GBND')

%exportgraphics(gca,'../pic/gbndDynamicUniformSim.png')

%%

plot(conv(inv(odfA),bndA),'complete','upper')
setColorRange([0,4.2])
exportgraphics(gca,'../pic/gbndDynamicUniform.png')

%% -------------------------------------------------------
%% Example 2 - textured microstructure

%f = fibre(Miller(1,1,1,cs),zvector)
%odfA = fibreODF(f)

odfA = fibreODF(fibre.gamma(cs))

plotPDF(odfA,Miller({1,0,0},cs))
setColorRange([0,3.3])
%mtexColorbar

%%

%modf = unimodalODF(moriRef,'halfwidth',1.5*degree)

odfB = rotate(odfA,moriRef,'right');
%odfB = conv(odfA,modf)

odf = conv(odfB + odfA,SO3DeLaValleePoussinKernel('halfwidth',2*degree)) /2;

plotPDF(odf,Miller({1,0,0},cs))
setColorRange([0,2.5])

%exportgraphics(gca,'../pic/pfDynamicTex.png')

%%

plotPDF(odfA,Miller({1,0,0},cs))
setColorRange([0 3.5])
%exportgraphics(gca,'../pic/odfA.png')

%%

plotPDF(odfB,Miller({1,0,0},cs))
setColorRange([0 3.5])
mtexColorbar
%exportgraphics(gca,'../pic/odfB.png')

%%
plotPDF(odf,Miller({1,0,0},cs))
%exportgraphics(gca,'../pic/odfB.png')



%% simulate the microstructure

neper.morpho = 'diameq:normal(10,1),aspratio(1,1,1)';
grains = neper.simulateTwinGrains(5000,moriRef,odf,bndA)

%save("grains32b.mat",'grains','moriRef','bndA','cs')

%%

load("grains32b.mat")

%%

plot(grains,grains.meanOrientation)

pC = plottingConvention.default3D;
pC.setView

%exportgraphics(gca,'../pic/grainDynamicTex.png')

%% measured ODF

odfABSim = calcDensity(grains.meanOrientation,'halfwidth',2.5*degree)

plotPDF(odfABSim,Miller({1,0,0},cs))
%setColorRange([0,2.2])
mtexColorbar

%exportgraphics(gca,'../pic/pfDynamicTexSim.png')

%% simulated GBCD
% since mori = inv(mori) we can not distinguish between gbndA and gbndB

gB = grains.boundary('indexed');

gbndABSim = gB.calcGBND(grains,moriRef,'halfwidth',1.5*degree)

figure(1)
plot(gbndABSim,'complete','upper')
mtexColorbar
mtexTitle('simulated crystalographic GBND')

%exportgraphics(gca,'../pic/gbndCrystDynamic.png')

%% simulated GBND

bndABsim = calcGBND(gB,grains,'halfwidth',1*degree)

plot(bndABsim,'complete','upper')


%% specimen GBND - simulated

gB = grains.boundary('indexed');
gbndS = gB.calcGBND('halfwidth',5*degree)

% simulated specimen GBND
plot(gbndS,'figSize','normal')
setColorRange([0,4.2])
%mtexTitle('experimental specimen GBND')
%exportgraphics(gca,'../pic/gbndDynamicTexSim.png')

%% predicted specimen GBND with theoretic odf and GBCD

plot(conv(inv(odfA),bndA))
setColorRange([0,4.2])
%mtexTitle('theoretic specimen GBND')
%setColorRange('equal')
%mtexColorbar

%exportgraphics(gca,'../pic/gbndDynamicTex.png')

%%

%nextAxis
plot(conv(inv(odfA),bndA))
nextAxis
plot(conv(inv(odfB),bndB))
%mtexTitle('theoretic specimen GBND')
%setColorRange('equal')
mtexColorbar

%%


% theoretic specimen GBND with simulated odf and GBCD
plot(conv(inv(odfABSim),bndABsim))
%mtexTitle('theoretic specimen GBND')
%setColorRange('equal')

mtexColorbar
%exportgraphics(gca,'../pic/gbndDynamicTexAB.png')

% saveFigure('../pic/gbndSpecimenDynamic.png')


%% -----------------------------------------------------------

odfSim = calcDensity(grains.meanOrientation,'weights',grains.volume,'halfwidth',10*degree)

plotPDF(odfSim,Miller(1,0,0,cs))

% exportgraphics(gca,'../pic/pdfUniformDynamic.png')

%% if it would be a macroscopically driven microstructure and following image would agree

nextAxis

plot(conv(odf,gbndS),'complete','upper')
mtexTitle('theoretic crystallographic GBND')
setColorRange('equal')