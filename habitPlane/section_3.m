% MTEX script for the paper
%
%% Habit plane determination from reconstructed parent phase orientation maps
%
% by Frank Niessen, Tuomo Nyyssönen, Azdiar A. Gazder, Ralf Hielscher:
%
% * to run the script at least MTEX 5.9 has to be installed
% 
%% Section 3.2 - Trace Determination
% 
%% Data Import 

load('mat/ebsd2.mat')

%% Parent Grain reconstruction
% In this section we focus on a single parent grain. For that reason lets
% first reconstruct the parent grain structure.

ebsdParent = job.calcParentEBSD;
ebsdParent = ebsdParent(job.csParent);
[parentGrains,ebsdParent('indexed').grainId] = calcGrains(ebsdParent('indexed'),'angle',5*degree,'qhull');
ebsdParent = ebsdParent.gridify;

%% The test parent grain colorized by variantId
% Next we select the largest parent grain and plot it together with all its
% child variants.

% select largest parent grain
[~,pId] = max(parentGrains.grainSize);

plot(ebsdParent(ebsdParent.grainId==pId),ind2color(ebsdParent(ebsdParent.grainId==pId).variantId))

[~,mId] = max(job.parentGrains.grainSize);
mId = job.parentGrains.id(mId);

hold on
plot(job.grainsPrior(job.mergeId ==  mId).boundary,'linewidth',1)
hold off

%% Figure 3 Fourier Transform based Trace Determination

% Fig. 3 (a)
figure(1) 
plot(parentGrains(pId).boundary)
hold on

% the variants we are interested in
vIds = [3,13];

for cId = 1:2

  vId = vIds(cId);

  % now we are performing the Fourier based trace determination

  % Step 1: black and white image of the positions of the variant
  X = ebsdParent.variantId==vId & ebsdParent.grainId == pId;
  ebsdLocal = ebsdParent.subGrid(X);
  X = trimMat(X);
  
  % Step 2:  Fourier transform
  Y = rescale(abs(fftshift(fft2(X - mean(X,'all')))));
  
  % cut of circle
  mxy = (1 + size(X))./2;
  dr = min(mxy)/2;
  [x,y] = meshgrid(1:size(Y,2),1:size(Y,1));
  
  % cut of circle
  r = sqrt((x-mxy(2)).^2+(y-mxy(1)).^2);
  Y = Y .* (r<dr);
  
  % Step 3: PCA
  M = Y(:) .* [x(:)-mxy(2),y(:)-mxy(1)];
  [lambda,v,~] = eig2(M' * M);
  t = vector3d(v(1),v(2),0,'antipodal');

  % Figure 3 (b / c)
  % display Fourier transform
  figure(1+cId)
  imagesc(x(1,:)-mxy(2),y(:,1)-mxy(1),flipud(sqrt(Y))), axis equal tight;
  hold on
  plot(-dr*v(1)*[-1,1],dr*v(2)*[-1,1],'LineWidth',2,'Color',ind2color(cId))
  hold off
  mtexColorMap white2black
  clim([0,0.55])
  xlim([-dr,dr])
  ylim([-dr,dr])
  axis off

  % Figure 3 (a)
  % add trace to the variant figure
  figure(1)
  plot(ebsdLocal,'facecolor',ind2color(cId),'micronbar','off')
  hold on
  quiver(ebsdLocal(ceil(end/2),ceil(end/2)),10*t,'color',ind2color(cId))
  legend off

end
hold off

%exportgraphics(gca,'../images/grain1V13Eig.png','Resolution',600)
%saveFigure('../images/grain1V15.png');
%exportgraphics(gca,'../images/grain1V1.png','Resolution',600)
%exportgraphics(gca,'../images/grain1V3_13.png','Resolution',600)
%exportgraphics(gca,'../images/grain1V3_13_g.png','Resolution',600)

%% Figure 4: Radon Transform Based Trace Detection

% Here we consider only one variant
vId = 3;

% Step 1: black and white image of the positions of the variant
X = ebsdParent.variantId==vId & ebsdParent.grainId == pId;
ebsdLocal = ebsdParent.subGrid(X);

% Step 2: Radon transform
R = radon(trimMat(X));

%% Fig. 4 (a)
% display the radon transform

imagesc(sqrt(rescale(R)))
mtexColorMap white2black
xlabel('angle in degree')
ylabel('shift n')
xlim([0,180])
xticks(0:30:180)
yticks([])

%% Fig 4 (b)

% Step 3: Fourier transform of the sinogram
Y = fftshift(abs(fft(R,[],1)),1);

imagesc(Y(1:floor(end/2)-3,:))
mtexColorMap white2black
xlabel('angle in degree')
ylabel('frequency')
xlim([0,180])
xticks(0:30:180)
yticks([])
% 
%axis off, exportgraphics(gca,'../images/grain1V3RadonF.png','Resolution',600)

%% Fig 4 (c)

% Step 4: sum over all frequencies
Z = sum(Y(1:floor(end/2)-3,:));

plot(1:180,Z,'linewidth',2);
xlabel('angle in degree')
xlim([0,180])
xticks(0:30:180)
yticks([])
%xticks([]), exportgraphics(gca,'../images/grain1V3RadonFSum.png','Resolution',600)

%% Fig. 5 (a) Boundary Segment based trace dection
% Illustration of grain boundary smoothing 

childGrains = job.grainsPrior

% select a region
reg = [232,238,47,53];
plot(ebsd('Fe-BCC'),ebsd('Fe-BCC').orientations,...
  'region',reg,'micronbar','off','FaceAlpha',0.5)

% boundaries before smoothing
hold on
plot(childGrains.boundary,'region',reg,'linewidth',2)
hold off

%exportgraphics(gca,'../images/grainSmoothing0.pdf')

%% Fig. 5 (b) circular histogram of boundary dirctions

% Fig 3(b)
N = 48;
h = histogram(job.grainsPrior.boundary.direction,'normalization','pdf',...
  'BinEdges',pi/N + linspace(0,2*pi,N+1));
ax = gca;
ax.ThetaTickLabel = [];
ax.RTickLabel = [];
%exportgraphics(gca,'../images/hist0.pdf')

%% Fig 5 (c) 

% smooth grain boundaries
gg = smooth(childGrains,10)

plot(ebsd('Fe-BCC'),ebsd('Fe-BCC').orientations,...
  'region',reg,'micronbar','off','FaceAlpha',0.5)

hold on
plot(gg.boundary,'region',reg,'linewidth',2)
hold off
%exportgraphics(gca,'../images/grainSmoothing10.pdf')

%% Fig 5(d) circular histogram of smoothed boundary dirctions

histogram(gg.boundary.direction,'BinEdges',pi/N + linspace(0,2*pi,N+1))
ax = gca;
ax.ThetaTickLabel = [];
ax.RTickLabel = [];
%exportgraphics(gca,'../images/hist10.pdf')

%% Fig 6 

% Fig 6 (a)
close all, figure(1)
plot(parentGrains(pId).boundary)
hold on

for cId = 1:2

  vId = vIds(cId);

  vGrains = gg(job.mergeId == mId & job.variantId == vId);
  gB = vGrains.boundary;
  
  % Fig 6 (a)
  figure(1), hold on
  plot(vGrains,'FaceColor',ind2color(cId),'FaceAlpha',0.2);

  % Figure 6 (b) - circular histogram
  figure(2)
  N = 64;
  histogram(gB.direction,'weights',gB.segLength,...
    'BinEdges',0+pi/N + linspace(0,2*pi,N+1),'facecolor',ind2color(cId));
  set(gca,'RTickLabel',[],'ThetaTickLabel', [])
  hold on
  
  % density function
  fun = calcDensity(gB.direction.rho, 'weights',gB.segLength,'periodic','sigma',10*degree);
  fun.antipodal = true;
  fun.fhat = 1.2*fun.fhat / N;
  plot(fun,'linewidth',3,'color',ind2color(cId));
  
  % Figure 6 (c) - characteristic shape
  figure(3)
  cS = characteristicShape(vGrains.boundary ,'noSimplify');
  [omega,a,b] = principalComponents(cS);
  t = vector3d.byPolar(pi/2,omega,'antipodal');
  
  plot(cS,'Color',ind2color(cId),'linewidth',2)
  set(gca,'RTickLabel',[],'ThetaTickLabel', [])
  hold on

  % Fig 6 (a)
  figure(1)
  hold on
  quiver(ebsdLocal(ceil(end/2),ceil(end/2)),10*t,'color',ind2color(cId))

end

hold off
legend off

%exportgraphics(gca,'../images/grain1V3_13_hist2.pdf','Resolution',600)
%exportgraphics(gca,'../images/grain1V3_13_shape_S.png','Resolution',600)
%saveFigure('../images/grain1V15.png');
