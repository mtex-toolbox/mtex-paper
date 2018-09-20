% Script file for the publication 
%
%% Denoising of Crystal Orientation Maps
% 
% by Ralf Hielscher, Christian B. Silbermann, and Eric Schmidl
%
%% Section 3.1  Kernel Average Misorientation Estimation
%

% where to save the image files
pname = ['..' filesep 'pic' filesep 'sim' filesep];
plotx2east


%% The KAM of the noise free data - Fig. 3.1 (left)

ebsd = simEBSD;
kam = ebsd.KAM('threshold',10*degree,'order',3)./degree;

plot(ebsd,kam,'micronbar','off')
CLim(gcm,[0,0.5])
mtexColorbar
%saveFigure([pname '/pic/kamSimNoiseFree.png'])



%% next we increase the noise and look for the KAM

d = kam(12,:).';
for hw = [0.02 0.1 0.2 0.5 1]
  
  ebsd = simEBSD('poussin',hw*degree);

  kam = ebsd.KAM('threshold',10*degree,'order',3)./degree;
  
  d(:,end+1) = kam(12,:);
  
end

d = [(1:100).',d];

%cprintf(d,'-fc',['../pic/KAM/kamSim.txt'],'-q',true);


%%  The KAM of noisy data - Fig. 3.1 (right)

ebsd = simEBSD('poussin',1*degree);
plot(ebsd,kam(:),'micronbar','off','figSize','tiny')
CLim(gcm,[0,4])
mtexColorbar
%saveFigure(['../pic/kamSimNoise.png'])


%% KAM profiles for different noise - Fig. 3.2

clf
plot(d(:,1),d(:,2:end),'linewidth',2)


%% KAM maps for different regularization parameter alpha - Fig. 3.4

F = halfQuadraticFilter;

d = (1:100).';
for a = [0,0.01,0.04,0.1,0.4,1]

  F.alpha = a;
  
  ebsd_smoothed = smooth(ebsd,F);

  kam = ebsd_smoothed.KAM('threshold',5*degree,'order',3)./degree;
  
  %plot(ebsd,kam(:))
  %saveFigure([pname '/pic/kamSim' xnum2str(1000*a) '.png'])
  
  %d(:,end+1) = nanmedian(kam);
  d(:,end+1) = kam(12,:);
  
end

%cprintf(d,'-fc',['../pic/KAM/kamSimSec.txt'],'-q',true);

%% KAM profiles for different regularization parameter alpha - Fig. 3.3

clf
plot(d(:,1),d(:,2:end))
