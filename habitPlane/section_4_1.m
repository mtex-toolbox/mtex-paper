% MTEX script for the paper
%
%% Habit plane determination from reconstructed parent phase orientation maps
%
% by Frank Niessen, Tuomo Nyyssönen, Azdiar A. Gazder, Ralf Hielscher:
%
% * to run the script at least MTEX 5.9 has to be installed
% 
%% Section  4.1 Sensitivity to deviating trace directions
%
%% Preperation

% define phases
cs_fcc = crystalSymmetry('432','mineral', 'fcc'); % parent
cs_bcc = crystalSymmetry('432','mineral', 'bcc'); % child

% Determine OR
p2c = orientation.GreningerTrojano(cs_fcc,cs_bcc);

% assumed habit plane in parent coordinates
h_fcc = Miller({5,7,5},cs_fcc,'antipodal');

%% Case 1: symmetry operators with no chance of error, traces have errors

% number of parent orientations 
l = 150;

% some random parent orientations
ori_fcc = orientation.rand(l,cs_fcc);

% and consider all variants
ori_fcc_var = ori_fcc.symmetrise;

% Get perfect trace vectors on imaging plane:
traceImPlane = cross(ori_fcc_var .* h_fcc,zvector);


%% Figure 7 (a) - simulate noisy traces on the imaging plane

% generate rotational angles by a Gaussian distribution
sigma = 12; % standard deviation
mu1 = 0;    % center
mu2 = 10;
omega1 = sigma.*randn(size(traceImPlane)) + mu1; % rotations on normal distribution
omega2 = sigma.*randn(size(traceImPlane)) + mu2; % rotations on normal distribution

figure
histogram(omega1,'Normalization','pdf')
hold on
histogram(omega2,'Normalization','pdf')
fplot(@(x) ...
    1/(sigma*sqrt(2*pi))*exp(-(x-mu1)^2/(2*sigma^2)),'linewidth',2,'color','black')

fplot(@(x) ...
    1/(sigma*sqrt(2*pi))*exp(-(x-mu2)^2/(2*sigma^2)),'linewidth',2,'color','black')
xlabel('Rotation angle, [°]')
legend(['\sigma = 12' newline '\mu = 0'],...
    ['\sigma = 12' newline '\mu = 10'],...
    'Location','northwest')

% Perform rotations to simulate noise:
traceImPlane1 = rotation.byAxisAngle(zvector, omega1*degree) .* traceImPlane;
traceImPlane2 = rotation.byAxisAngle(zvector, omega2*degree) .* traceImPlane;

%% Recover habit plane from the noisy traces

% Rotate noisy traces back to the parent orientation reference frame
traceParRefFrame1 = inv(ori_fcc_var) .* traceImPlane1;
traceParRefFrame2 = inv(ori_fcc_var) .* traceImPlane2;

% Determine the habit plane from the rotated traces
hPlane1 = perp(traceParRefFrame1,'robust');
hPlane2 = perp(traceParRefFrame2,'robust');

round(hPlane1)
round(hPlane2)

%% Section 4.2: Sensitivity to misindexed symmetry operations

% compute perfect child orientations
ori_bcc = ori_fcc_var .* inv(p2c);

% compute variant and packet id
[vId,pId] = calcVariantId(ori_fcc_var,ori_bcc,p2c);

% twinning planes axes ordered by packet Id
h1 = Miller({1,1,1},{1,-1,1},{-1,1,1},{1,1,-1},cs_fcc,'antipodal');

% generate twinned parent orientations that share most of the variants
twinning = orientation.byAxisAngle(h1(pId),60*degree);

% apply twinning to 10 percent of all parent orientations
ind = rand(size(ori_bcc))>0.75;
ori_fcc_twin = ori_fcc_var;
ori_fcc_twin(ind) = ori_fcc_twin(ind).project2FundamentalRegion .* twinning(ind).';

% compute the parentId from the wrong parent orientations
vIdWrong = calcVariantId(ori_fcc_twin, ori_bcc, p2c);
vIdWrong = reshape(vIdWrong,size(ori_bcc));

% and use it for projecting the traces into parent coordinates
traceParRefFrame3 = inv(p2c) .* variants(p2c,vIdWrong) ...
    .* inv(ori_fcc_twin.project2FundamentalRegion) .* traceImPlane;

hPlane3 = perp(traceParRefFrame3,'robust');

round(hPlane3)

scatter(traceParRefFrame3,'MarkerFaceColor','none')

%% Section 4.2 (ii): habit planes with a poor representative OR

% add some noise to the child orientations
ori_bcc = ori_fcc_var .* inv(p2c);
ori_bcc_rotated = rotation.rand(size(ori_bcc),'maxAngle',10*degree) .* ori_bcc;

% compute variantId with respect to a non representative OR - here KS
p2c_KS = orientation.KurdjumovSachs(cs_fcc,cs_bcc);

[vId_KS] = calcVariantId(ori_fcc_var,ori_bcc_rotated,p2c_KS);
vId_KS = reshape(vId_KS,size(ori_bcc));

traceParRefFrame4 = inv(p2c_KS) .* variants(p2c_KS,vId_KS) ...
    .* inv(ori_fcc_var.project2FundamentalRegion) .* traceImPlane;

hPlane4 = perp(traceParRefFrame4,'robust');

round(hPlane4)

%% Case 4: habit planes with a poor representative OR and traces with errors

traceParRefFrame5 = symrots_true(vId_KS) ...
    .*inv(ori_fcc(z).project2FundamentalRegion).*traceImPlane1;

hPlane5 = perp(traceParRefFrame5,'robust');

round(hPlane5)


%% Quantiles and median angles between habit and traces
delta = 90 - angle(hPlane1,traceParRefFrame1,'noSymmetry')/degree;
% quant = mean(delta,[0.25 0.5 0.75])
quant = mean(delta)
delta2 = 90 - angle(hPlane2,traceParRefFrame2,'noSymmetry')/degree;
% quant2 = quantile(delta2,[0.25 0.5 0.75])
quant2 = mean(delta2)
delta3 = 90 - angle(hPlane3,traceParRefFrame3,'noSymmetry')/degree;
% quant3 = quantile(delta3,[0.25 0.5 0.75])
quant3 = mean(delta3)
delta4 = 90 - angle(hPlane4,traceParRefFrame4,'noSymmetry')/degree;
% quant4 = quantile(delta4,[0.25 0.5 0.75])
quant4 = mean(delta4)
delta5 = 90 - angle(hPlane5,traceParRefFrame5,'noSymmetry')/degree;
% quant5 = quantile(delta5,[0.25 0.5 0.75])
quant5 = mean(delta5)
%% Angle between original plane and reconstructed habit plane

ang1 = angle(hPlane1,h_fcc)/degree
ang2 = angle(hPlane2,h_fcc)/degree
ang3 = angle(hPlane3,h_fcc)/degree
ang4 = angle(hPlane4,h_fcc)/degree
ang5 = angle(hPlane5,h_fcc)/degree

%%
close all

setMTEXpref('FontSize',20);

plotx2east
plotzIntoPlane

% Plot the traces and the fitted habit plane:
figure
mtexFigure('figsize','medium')
scatter(traceParRefFrame1,'MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.1,'projection','eangle')
hold on
circle(hPlane1,'linewidth',6,'linecolor','w')
circle(hPlane1,'linewidth',4,'linecolor','k')
% circle(h_fcc,'linewidth',3,'linecolor','r')

% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame1,'noSymmetry'),'projection','eangle')
% circle(hPlane1,'linewidth',6,'linecolor','w')

figure
mtexFigure('figsize','medium')
scatter(traceParRefFrame2,'MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.1,'projection','eangle')
hold on
circle(hPlane2,'linewidth',6,'linecolor','w')
circle(hPlane2,'linewidth',4,'linecolor','k')

% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame2,'noSymmetry'),'projection','eangle')
% circle(hPlane2,'linewidth',6,'linecolor','w')
%%
figure
mtexFigure('figsize','medium')
plot(hPlane3,'linewidth',10,'linecolor','r','projection','eangle','circle')
    circle(Miller(1,0,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(-1,0,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(0,-1,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(0,1,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(1,1,0,cs_fcc),'color','k','linewidth',2)
    circle(Miller(1,-1,0,cs_fcc),'color','k','linewidth',2)
    circle(Miller(1,0,0,cs_fcc),'color','k','linewidth',2)
    circle(Miller(0,1,0,cs_fcc),'color','k','linewidth',2)
hold on
scatter(traceParRefFrame3,'MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.1,'projection','eangle', ...
    'MarkerFaceColor',[0 0.4470 0.7410])
plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
    'MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',4)
plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
    'MarkerFaceColor','none','MarkerEdgeColor','w','linewidth',2)


% %%
% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame3,'noSymmetry'),'projection','eangle')
% circle(hPlane3,'linewidth',6,'linecolor','w')
% circle(symrots_true*(h_fcc),'linewidth',2,'linecolor','k')

%%
figure
mtexFigure('figsize','medium')
plot(hPlane4,'linewidth',10,'linecolor','r','projection','eangle','circle',...
    'LineEdgeAlpha')
hold on
scatter(traceParRefFrame4,'MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.1,'projection','eangle', ...
    'MarkerFaceColor',[0 0.4470 0.7410])
hold on

%     circle(symmetrise(Miller(1,1,1,cs_fcc)),'color','k','linewidth',2)
    circle(Miller(1,0,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(-1,0,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(0,-1,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(0,1,1,cs_fcc),'color','k','linewidth',2)
    circle(Miller(1,1,0,cs_fcc),'color','k','linewidth',2)
    circle(Miller(1,-1,0,cs_fcc),'color','k','linewidth',2)
    circle(Miller(1,0,0,cs_fcc),'color','k','linewidth',2)
    circle(Miller(0,1,0,cs_fcc),'color','k','linewidth',2)
hold on
plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
    'MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',4)
plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
    'MarkerFaceColor','none','MarkerEdgeColor','w','linewidth',2)
% 
%     circle(Miller(5,7,5,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(7,5,5,cs_fcc),'color','k','linewidth',2)

%
% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame4,'noSymmetry','halfwidth',5*degree))
% hold on
% plot(hPlane3,'linewidth',4,'linecolor','r','projection','eangle','circle',...
%     'LineEdgeAlpha')
% %     circle(symmetrise(Miller(1,1,1,cs_fcc)),'color','k','linewidth',2)
%     circle(Miller(1,0,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(-1,0,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(0,-1,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(0,1,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(1,1,0,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(1,-1,0,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(1,0,0,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(0,1,0,cs_fcc),'color','k','linewidth',2)
% hold on
% plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
%     'MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',4)
% plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
%     'MarkerFaceColor','none','MarkerEdgeColor','w','linewidth',2)
% 
% text(hPlane4,['Quartiles:' newline ...
%     num2str(round(quant4(1),1)) ' °' newline ...
%     num2str(round(quant4(2),1)) ' °' newline ...
%     num2str(round(quant4(3),1)) ' °' newline ...
%     'angle = ' num2str(round(ang4,1)) ' °'],'color','k',...
%     'Interpreter','latex')

%%
% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame4,'noSymmetry'),'projection','eangle')
% circle(hPlane4,'linewidth',6,'linecolor','w')
% circle(symrots_true*(h_fcc),'linewidth',2,'linecolor','k')

%%
figure
mtexFigure('figsize','medium')
scatter(traceParRefFrame5,'MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.1,'projection','eangle')
hold on
circle(hPlane5,'linewidth',6,'linecolor','w')
circle(hPlane5,'linewidth',4,'linecolor','k')

%     circle(Miller(5,7,5,cs_fcc),'color','k','linewidth',2)
%      circle(Miller(7,5,5,cs_fcc),'color','k','linewidth',2) 

% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame5,'noSymmetry'))
% hold on
% plot(hPlane3,'linewidth',4,'linecolor','r','projection','eangle','circle',...
%     'LineEdgeAlpha')
% %     circle(symmetrise(Miller(1,1,1,cs_fcc)),'color','k','linewidth',2)
%     circle(Miller(1,0,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(-1,0,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(0,-1,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(0,1,1,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(1,1,0,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(1,-1,0,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(1,0,0,cs_fcc),'color','k','linewidth',2)
%     circle(Miller(0,1,0,cs_fcc),'color','k','linewidth',2)
% hold on
% plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
%     'MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',4)
% plot(Miller(-1,0,1,cs_fcc),'Marker','square','MarkerSize',15,...
%     'MarkerFaceColor','none','MarkerEdgeColor','w','linewidth',2)
%%
% figure
% mtexFigure('figsize','medium')
% plot(calcDensity(traceParRefFrame5,'noSymmetry'),'projection','eangle')
% circle(hPlane5,'linewidth',6,'linecolor','w')
% circle(symrots_true*(h_fcc),'linewidth',2,'linecolor','k')
