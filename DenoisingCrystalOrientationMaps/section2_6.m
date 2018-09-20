% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 2.6 Variational Filters
%

% where to save the image files
pname = ['..' filesep 'pic' filesep 'sim' filesep];
plotx2east

% de la Vallee Poussin distributed noise
ebsd = simEBSD('poussin',2*degree);

d(:,1) = 1:100;

%% --------------- smoothing splines ----------------------------------

F = splineFilter;

%% alpha = 0.2, Fig. 2.5a

F.alpha = 0.2;

ebsd_smoothed = smooth(ebsd,F)
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure([pname 'splineSim1.png'])
d(:,2) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% alpha = 5, Fig. 2.5b

F.alpha = 5;

ebsd_smoothed = smooth(ebsd,F);

plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')
%saveFigure([pname 'splineSim2.png'])
d(:,3) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% alpha = 5 and impulsive noise Fig. 2.5c

ebsd = simEBSD('poussin',2*degree,'saltPepper',0.05);

F.alpha = 5;
ebsd_smoothed = smooth(ebsd,F)

plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
caxis([0,18])
%saveFigure([pname 'splineSim3.png'])
d(:,4) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%cprintf(d,'-fc',[pname 'splineSimSec.txt'],'-q',true);

%% --------------------- total variation -----------------------------------

F = halfQuadraticFilter;
ebsd = simEBSD('poussin',2*degree);

%% alpha = 0.025, Fig. 2.6a

F.alpha = 0.025;

ebsd_smoothed = smooth(ebsd,F)
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure([pname 'splineSim1.png'])
d(:,2) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% alpha = 0.5, Fig. 2.6b

F.alpha = 0.5;

ebsd_smoothed = smooth(ebsd,F);

plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')
%saveFigure([pname 'splineSim2.png'])
d(:,3) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% alpha = 0.5 and impulsive noise Fig. 2.6c

ebsd = simEBSD('poussin',2*degree,'saltPepper',0.05);

ebsd_smoothed = smooth(ebsd,F)

plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
caxis([0,18])
%saveFigure([pname 'splineSim3.png'])
d(:,4) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%cprintf(d,'-fc',[pname 'splineSimSec.txt'],'-q',true);

%% --------------------- infimal convolution -----------------------------------

F = infimalConvolutionFilter;
ebsd = simEBSD('poussin',2*degree);

%% alpha = 0.025, Fig. 2.7a

F.lambda = 0.025;
%F.lambda = 0.1;
F.mu = 0.05;

ebsd_smoothed = gridify(smooth(ebsd,F));
plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure([pname 'splineSim1.png'])
d(:,2) = angle(ebsd_smoothed(12,:).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%plot(d(:,3))

%% alpha = 0.5, Fig. 2.7b

F.lambda = 0.1;
F.mu = 0.2;

ebsd_smoothed = gridify(smooth(ebsd,F));

plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')
%saveFigure([pname 'splineSim2.png'])
d(:,3) = angle(ebsd_smoothed(12,:).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% alpha = 0.5 and impulsive noise Fig. 2.7c

F.lambda = 0.025;
%F.lambda = 0.1;
F.mu = 0.05;

ebsd = simEBSD('poussin',2*degree,'saltPepper',0.05);

ebsd_smoothed = smooth(ebsd,F)

plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
caxis([0,18])
%saveFigure([pname 'splineSim3.png'])
d(:,4) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%cprintf(d,'-fc',[pname 'splineSimSec.txt'],'-q',true);

