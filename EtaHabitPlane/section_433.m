% Script file for section 4.3.3 of the publication 
%
%% Interplay between habit plane and orientation relationship 
%% using the example of eta'-Al8Fe3 in eta-Al5Fe 2:
%  an electron backscatter analysis
% 
% by Hanka Becker, Ralf Hielscher and Andreas Leineweber
%
% For running this script MTEX 5.8 is required
%
%% 1. determine parent and child variant ids
%
% p2c is the reference OR
candidates = (csChild * p2c * csParent).';

% child x parent;
d = angle_outer(candidates,p2cmori,'noSymmetry')./degree;

[~,id] = min(d);
[parentId,childId] = ind2sub(size(candidates),id);

%% Figure 9a - Child traces with variant id

scatter(traceChild,ind2color(childId))

vChild = Miller({1,0,0},{0,1,0},{0,0,1},csChild,'uvw');
annotate(vChild,'labeled','backgroundColor','w')

hold on
circle(nChild.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceChildVariant.png')

%% Figure 9d - Parent traces with variant id

scatter(traceParent,ind2color(parentId))

annotate(vParent,'labeled','backgroundColor','w')

hold on
circle(nParent.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceParentVariant.png')

%% reassign variants id

% pseudo symmetry correction
traceChildCor = reshape(csChild.rot(childId),[],1) .* traceChild;

% omega  - fit with the guessed normal
% omega2 - fit with the rotated normal
omega = pi/2-angle(nChild,traceChildCor,'noSymmetry');
omega2 = pi/2-angle(nChild,csChild.rot(2)*traceChildCor,'noSymmetry');
ind = omega > omega2(:);

% reassign variants ids
parentId(ind) = 5 - parentId(ind);
childId(ind) = 3 - childId(ind);

%% Figure 9b - Child traces with corrected variant id

scatter(traceChild,ind2color(childId))

vChild = Miller({1,0,0},{0,1,0},{0,0,1},csChild,'uvw');
annotate(vChild,'labeled','backgroundColor','w')

hold on
circle(nChild.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceChildVariantCor.png')

%% Figure 9e - Parent traces with corrected variant id

scatter(traceParent,ind2color(parentId))

annotate(vParent,'labeled','backgroundColor','w')

hold on
circle(nParent.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceParentVariantCor.png')

%% project traceParent and traceChild onto its fir

traceParentCor = reshape(csParent.rot(parentId),[],1) .* traceParent;
traceChildCor = reshape(csChild.rot(childId),[],1) .* traceChild;

%% Determine the best fitting plane normal in child coordinates

% a search grid of possible normal directions
n = plotS2Grid('resolution',0.125*degree,'upper');

chi = mean((pi/2-abs(angle_outer(traceChildCor,n,'noSymmetry'))).^p,1).^(1/p);

[m,id1] = min(chi(:));

nChild = Miller(n(id1),csChild,'antipodal');
disp(['h_child = ' char(nChild ./ nChild.h)])
disp(['chi = ' xnum2str(m./degree)])

%% Figure 9c) - the traces in child coordinates 

scatter(traceChildCor)
annotate(vChild,'labeled','backgroundColor','w')

hold on
circle(nChild,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceChildProj.png')


%% Figure 9f) - The functional chi^2 in child coordinates 

contourf(Miller(n,csChild),chi./degree,'contours',20)
mtexColorMap black2white
mtexColorbar
%setColorRange([min(d(:)),10])

hold on
plot(nChild,'MarkerColor',ind2color(3))
hold off

%saveFigure('../../pic/chiChildProj.png')

%% Determine the best fitting plane normal in parent coordinates

% a search grid of possible normal directions
n = plotS2Grid('resolution',0.125*degree,'upper');

chi = mean((pi/2-abs(angle_outer(traceParentCor,n,'noSymmetry'))).^p,1).^(1/p);

[m,id1] = min(chi(:));

nParent = Miller(n(id1),csChild,'antipodal');
disp(['h_child = ' char(nParent ./ nParent.h)])
disp(['chi = ' xnum2str(m./degree)])


%% Figure 9d) - the projected traces in parent coordinates 

scatter(traceParentCor)
annotate(vParent,'labeled','backgroundColor','w')

hold on
circle(nParent,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceParentProj.png')


%% Figure 9xx) - The functional chi^2 in parent coordinates 

pcolor(Miller(n,csParent),chi./degree,'contours',20)
mtexColorMap black2white
mtexColorbar
%setColorRange([min(d(:)),10])

hold on
plot(nParent,'MarkerColor',ind2color(3))
hold off

%saveFigure('../../pic/chiParentProj.png')

%% Figure 9xx) - The functional chi_01 in parent coordinates 

pcolor(Miller(n,csParent),chi./degree,'contours',20)
mtexColorMap black2white
mtexColorbar
%setColorRange([min(d(:)),10])

hold on
plot(nParent,'MarkerColor',ind2color(3))
hold off

%saveFigure('../../pic/chi01ParentProj.png')
