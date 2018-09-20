% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 2.1 An artificial Orientation Map
%

% where to save the image files
pname = ['..' filesep 'pic' filesep 'sim' filesep];
plotx2east

%% some artifical data without noise - Fig. 2.1a

ebsd = simEBSD;

figure(1)
plot(ebsd,angle(ebsd.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
CLim(gcm,[0,20])
% saveFigure([pname 'sim.png'])

%%

figure(2)
xline = 12;
plot(angle(ebsd(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector)./degree)

d(:,1) = 1:100;
d(:,2) = angle(ebsd(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;


%% de la Vallee Poussin distributed noise  - Fig. 2.1b

ebsd = simEBSD('poussin',2*degree);

figure(1)
plot(ebsd,angle(ebsd.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
CLim(gcm,[0,20])
% saveFigure([pname 'poussin.png'])

%%

figure(2)
xline = 12;

plot(angle(ebsd(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector)./degree)
d(:,3) = angle(ebsd(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

%% salt and pepper noise - Fig. 2.1c


ebsd = simEBSD('poussin',2*degree,'saltPepper',0.05);

figure(1)
plot(ebsd,angle(ebsd.orientations*Miller(0,0,1,cs),zvector)./degree,'micronbar','off')
CLim(gcm,[0,20])
% saveFigure([pname 'saltPepper.png'])

%%
figure(2)
xline = 12;

plot(angle(ebsd(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector)./degree)
d(:,4) = angle(ebsd(xline*100+(1:100)).orientations*Miller(0,0,1,cs),zvector) ./ degree;

% store profiles in a text file
%cprintf(d,'-fc',[pname 'raw.txt'],'-q',true);