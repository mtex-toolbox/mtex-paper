% MTEX script for the paper
%
%% Modeling and Testing Grain Boundary Plane Distributions
%
% by Ralf Hielscher, Rüdiger Kilian, Katharina Tinka Marquardt, Erik Wünsche
%
% * to run the script at least MTEX 6.1 has to be installed
% 
%% Section 2.4 - Determining Boundary Distribution Functions from 3D-Microstructures
% 

load("grains31.mat")

%%

omega = linspace(-50,50,500);

figure(1)

psi = S2DeLaValleePoussinKernel('halfwidth',2.5*degree,'DisplayName','2.5°');

plotWithGap(omega, psi.eval(cos(omega*degree)),[400 1200],'ShowGapLabel',false,'linewidth',3)

hold on
psi = S2DeLaValleePoussinKernel('halfwidth',5*degree,'DisplayName','5°');
ax = gca;
ax.ColorOrderIndex = 2;
plot(omega, psi.eval(cos(omega*degree)),'LineWidth',3)
psi = S2DeLaValleePoussinKernel('halfwidth',10*degree,'DisplayName','5°');
plot(omega, psi.eval(cos(omega*degree)),'LineWidth',3)
hold off
 %
% psi = S2DeLaValleePoussinKernel('halfwidth',7.5*degree,'DisplayName','10°');

%hold on
%plot(psi,'linewidth',2,'symmetric')
%hold off

%xlim([-50 50])

l = legend('hw=2.5°','hw=5°','hw=10°','FontSize',16)

%exportgraphics(gcf,'../pic/psi.pdf', 'ContentType','vector')
%saveFigure('../pic/psi.pdf')


%%

figure(2)
gbnd = calcGBND(grains.boundary('indexed'),'kernel',psi)
plot(gbnd)
setColorRange([0,10.5])

%exportgraphics(gca,'../pic/gbnd_25.png')
%exportgraphics(gca,'../pic/gbnd_5.png')
%exportgraphics(gca,'../pic/gbnd_10.png')

%%

load("grains32a.mat")

%%
gB = discreteSample(grains.boundary('indexed'),1000)

%%


gbndA = calcGBND(gB,grains,'halfwidth',2.5*degree)

plot(gbndA,'complete','upper')

%exportgraphics(gca,'../pic/gbndA.png')


%%


gbcdnoSym = calcGBND(gB,grains,moriRef,'halfwidth',2.5*degree,'noSymmetry')

figure(1)
plot(gbcdnoSym,'complete','upper')
setColorRange([0,110])
mtexColorbar

%exportgraphics(gca,'../pic/gbcdnoSym.png')


%%

gbcd = calcGBND(gB,grains,moriRef,'halfwidth',2.5*degree)
figure(2)
plot(gbcd,'complete','upper')
setColorRange([0,110])
mtexColorbar


% exportgraphics(gca,'../pic/gbcdSym.png')