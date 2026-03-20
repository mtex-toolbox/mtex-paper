%% setup

cs = crystalSymmetry('432');

% Sigma 3 twinning
moriRef = orientation.byAxisAngle(Miller(1,1,1,cs),180*degree);

% the habit theoretic plane distribution
bndA = 0.4 * S2FunRBF(Miller(1,1,1,cs),'halfwidth',2.5*degree) ...
  + 0.6*S2FunRBF(Miller(1,0,1,cs),'halfwidth',5*degree);
bndA.antipodal = true;
bndA.CS = cs;

% compute intertwining symmetry elements
[~,~,csRed] = project2FundamentalRegion(moriRef,moriRef,moriRef);

bndB = rotate(bndA,inv(moriRef));
bndAB = 0.5*(bndA + bndB);

plot(bndA.symmetrise(csRed),'upper')
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)

%% simulate microstructure

neper.init
neper.geometry = "cube(4,4,4)"
neper.morpho = 'diameq:normal(10,1),aspratio(8,1,1)';

%odfA = fibreODF(fibre.gamma(cs))
odfA = uniformODF(cs);

grains = neper.simulateTwinGrains(4000,moriRef,odfA,bndA,'width',0.5)

% save("grains4Uniform.mat",'grains','moriRef','bndA','cs')

%%

load("grains4Uniform.mat")

%%

plot(grains,grains.meanOrientation)

pC = plottingConvention.default3D;
pC.setView
% exportgraphics(gca,'../pic/grainMixedUniform.png')

%% ODF estimation

odfEst = calcDensity(grains.meanOrientation,'weights',grains.volume,'halfwidth',10*degree)

%%

plotPDF(odfEst,Miller(1,0,0,cs))
%setColorRange([0,3.3])
mtexColorbar

% exportgraphics(gca,'../pic/pfMixSimUniform.png')

%% measured specimen GBND

gB = grains.boundary('indexed');
gbndS = gB.calcGBND('halfwidth',10*degree)

plot(gbndS,'figSize','normal')
mtexColorbar
%setColorRange([0,4.2])
%mtexTitle('experimental specimen GBND')

%exportgraphics(gca,'../pic/gbndMixTexSim.png')

%% measured crystal GBND

gbcd = gB.calcGBND(grains,moriRef,'halfwidth',2.5*degree)

figure(1)
plot(gbcd,'complete','upper')
mtexColorbar
%setColorRange([0,130])
%mtexTitle('simulated crystalographic GBND')

% exportgraphics(gca,'../pic/gbndABMixedSimUniform.png')

%% measured GBCD

gbndAB = gB.calcGBND(grains,'halfwidth',5*degree)

figure(1)
plot(gbndAB,'complete','upper')
mtexColorbar
%setColorRange([0,130])
%mtexTitle('simulated crystalographic GBND')

% exportgraphics(gca,'../pic/gbndABMixedSimUniform.png')



%% theoretical specimen GBND

plot(conv(inv(odfEst),gbndAB),'complete','upper')
mtexColorbar
setColorRange([0.95,1.05])
%exportgraphics(gca,'../pic/gbndABMixedUniform.png')





%% theoretical crys

plot(conv(inv(odfEst),gbndAB),'complete','upper')
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)
%mtexColorbar
setColorRange([0,130])

% exportgraphics(gca,'../pic/gbcdMixedUniform.png')

%% textured microstructure
% -------------------------------------------------------

load("grains4Tex.mat")

%%

plot(grains,grains.meanOrientation)

pC = plottingConvention.default3D;
pC.setView

%exportgraphics(gca,'../pic/grainMixedTex.png')

%% ODF estimation

odfEst = calcDensity(grains.meanOrientation,'weights',grains.volume,'halfwidth',5*degree)

%%

plotPDF(odfEst,Miller(1,0,0,cs))
%setColorRange([0,2.3])
mtexColorbar

% exportgraphics(gca,'../pic/pfMixTexSim.png')

%% simulated specimen GBND

gB = grains.boundary('indexed');
gbndS = gB.calcGBND('halfwidth',10*degree)

plot(gbndS,'figSize','normal')
mtexColorbar
%setColorRange([0,4.2])
%mtexTitle('experimental specimen GBND')

%exportgraphics(gca,'../pic/gbndMixTexSim.png')

%% theoretical specimen GBND

plot(conv(inv(odfEst),bndA),'complete','upper')
mtexColorbar
%setColorRange([0,4.2])
%exportgraphics(gca,'../pic/gbndMixedTex.png')


%% simulated crystal GBCD

gB = grains.boundary('indexed');

gbndAB = gB.calcGBND(grains,moriRef,'halfwidth',2.5*degree)

figure(1)
plot(gbndAB,'complete','upper')
mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)
mtexColorbar
setColorRange([0,130])
%mtexTitle('simulated crystalographic GBND')

% exportgraphics(gca,'../pic/gbcdMixedSimUniform.png')

%% simulated crystal GBND

gbndAB = gB.calcGBND(grains,'halfwidth',5*degree)

plot(gbndAB,'complete','upper')
%mtexTitle('$60^{\circ}/\left<111\right>$','fontSize',20)
mtexColorbar
%setColorRange([0,130])

% exportgraphics(gca,'../pic/gbndABMixedTexSim.png')


%% theoretical crystal GBND

plot(conv(odfEst,gbndS),'complete','upper')
mtexColorbar
%setColorRange([0,130])

% exportgraphics(gca,'../pic/gbndABMixedTex.png')

%%






%% simulated GBCD
% since mori = inv(mori) we can not distinguish between gbndA and gbndB

gB = grains.boundary('indexed');

gbndAB = gB.calcGBND(grains,moriRef,'halfwidth',2.5*degree)

figure(1)
plot(gbndAB,'complete','upper')
mtexColorbar
mtexTitle('simulated crystalographic GBND')

%exportgraphics(gca,'../pic/gbndCrystDynamic.png')

%% simulated GBND

bndABsim = calcGBND(gB,grains,'halfwidth',1.25*degree)

plot(bndABsim,'complete','upper')


%% specimen GBND

gB = grains.boundary('indexed');
gbndS = gB.calcGBND('halfwidth',5*degree)

% simulated specimen GBND
plot(gbndS,'figSize','normal')
setColorRange([0,4.2])
%mtexTitle('experimental specimen GBND')
%exportgraphics(gca,'../pic/gbndDynamicTexSim.png')

%%

% theoretic specimen GBND with theoretic odf and GBCD
plot(conv(inv(odf),bndA))
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
plot(conv(inv(odfABsim),bndABsim))
mtexTitle('theoretic specimen GBND')
%setColorRange('equal')
mtexColorbar

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


