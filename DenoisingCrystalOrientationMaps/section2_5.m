% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 2.5 The Kuwahara Filter
%

% where to save the image files
pname = ['..' filesep 'pic' filesep 'sim' filesep];
plotx2east

% de la Vallee Poussin distributed noise
ebsd = simEBSD('poussin',2*degree);

F = KuwaharaFilter;
d(:,1) = 1:100;

%% firster order neighbours Fig. 2.4a

F.numNeighbours = 1;

ebsd_smoothed = smooth(ebsd,F)
figure, plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')

%saveFigure([pname 'kuwaSim1.png'])
d(:,2) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% third order neighbours Fig. 2.4b

F.numNeighbours = 3;

ebsd_smoothed = smooth(ebsd,F);

figure, plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector),'micronbar','off')
%saveFigure([pname 'kuwaSim2.png'])
d(:,3) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% third order neighbours and impulsive noise Fig. 2.4c

ebsd = simEBSD('poussin',2*degree,'saltPepper',0.05);

ebsd_smoothed = smooth(ebsd,F)
figure, plot(ebsd_smoothed,angle(ebsd_smoothed.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
caxis([0,18])
%saveFigure([pname 'kuwaSim3.png'])
d(:,4) = angle(ebsd_smoothed(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%cprintf(d,'-fc',[pname 'kuwaSimSec.txt'],'-q',true);

