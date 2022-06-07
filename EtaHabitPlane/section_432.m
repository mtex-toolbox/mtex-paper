% script file for the publication 
%
%% Interplay between habit plane and orientation relationship 
%% using the example of eta'-Al8Fe3 in eta-Al5Fe 2:
%  an electron backscatter analysis
% 
% by Hanka Becker, Ralf Hielscher and Andreas Leineweber
%
% For running this script MTEX 5.8 is required
%
%% Section 4.3.2

%% data import

% define parent and child symmetry
csParent = crystalSymmetry('222', [7.6596 6.407 4.2344], 'mineral',  ...
'Al5_6Fe2_Cmcm_250_36', 'color', 'light blue');

csChild = crystalSymmetry('121', [7.6606 6.4244 4.1633],  ...
  [90,90.46,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral',  ...
  'Al5_6Fe2_C2c_250_365', 'color', 'light green');

% import orientation and trace data
chi = importdata('Trace-angles_Euler_new.xlsx');

chi = chi.data;
try chi = chi.Tabelle1; end

% parent and child orientations
oriChild = orientation.byEuler(chi(:,5:7).*degree,csChild);
oriParent = orientation.byEuler(chi(:,2:4).*degree,csParent);

% traces in specimen coordinates
trace = vector3d.byPolar(pi/2,chi(:,1)*degree,'antipodal');

% the traces with respect to the parent and the child
traceParent = inv(oriParent) .* trace;
traceChild = inv(oriChild) .* trace;

%% Section 4.2.2: Determination of the Orientation Relationship

p2cmori = inv(oriChild) .* oriParent;
p2c = project2FundamentalRegion(mean(p2cmori,'robust'))

matrix(p2c)

%% Section 4.3.2

%% Determine the best fitting plane normal in parent coordinates

% a search grid of possible plane normals
n = plotS2Grid('resolution',0.25*degree,'upper');

p = 2; % p = 0.1;

% compute chi
chi = mean((pi/2-abs(angle_outer(traceParent,n,'min'))).^p,1).^(1/p);

% find minumum value chi
[m,id1] = min(chi(:));

nParent = Miller(n(id1),csParent,'antipodal');
disp(['h_parent = ' char(nParent ./ nParent.h)])
disp(['chi = ' xnum2str(m./degree)])

%% Figure 8a) - the traces in parent coordinates

figure(1)
scatter(traceParent)
%scatter(traceParent,'label',1:length(traceParent))

vParent = Miller({1,0,0},{0,1,0},{0,0,1},csParent,'uvw');
annotate(vParent,'labeled','backgroundColor','w')

hold on
circle(nParent.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceParent.png')

%% Figure 8d) - The functional chi^2 in parent coordinates

figure(2)
contourf(Miller(n,csParent),chi ./ degree,'contours',20)
mtexColorMap black2white
mtexColorbar

hold on
plot(nParent.symmetrise,'MarkerColor',ind2color(3))
hold off

%saveFigure('../../pic/chiParent.png')

%% Determine the best fitting plane normal in child coordinates

% a search grid of possible normal directions
n = plotS2Grid('resolution',0.125*degree,'upper');

p = 2;
chi = mean((pi/2-abs(angle_outer(traceChild,n,'min'))).^p,1).^(1/p);

[m,id1] = min(chi(:));

nChild = Miller(n(id1),csChild,'antipodal');
disp(['h_child = ' char(nChild ./ nChild.h)])
disp(['chi = ' xnum2str(m./degree)])

%% Figure 8b) - the traces in child coordinates / uncorrected

scatter(traceChild)

vChild = Miller({1,0,0},{0,1,0},{0,0,1},csChild,'uvw');
annotate(vChild,'labeled','backgroundColor','w')

hold on
circle(nChild.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceChild.png')

%% Figure 8e) - The functional chi^2 in child coordinates / uncorrected

contourf(Miller(n,csChild),chi ./ degree,'contours',20)
mtexColorMap black2white
mtexColorbar
%setColorRange([min(d(:)),10])

hold on
plot(nChild.symmetrise,'MarkerColor',ind2color(3))
hold off

%saveFigure('../../pic/chiChild.png')

%% fit between parent and child habit plane normal 

min(angle(p2c * nParent.symmetrise, nChild)./degree)


%% correct child orientations

omega = pi/2-angle(nChild,traceChild,'min');
omega2 = pi/2-angle(nChild,csParent.rot(2).*traceChild,'min');

scatter(traceChild,omega./degree - omega2./degree,'label',1:43)

cond = omega2+5*degree < omega;
oriChild(cond) = oriChild(cond) * csParent.rot(2);

disp([' number fo orientations corrected: ', xnum2str(sum(cond))])

%% update child traces and OR p2c

traceChild = inv(oriChild) .* trace;
p2cmori = inv(oriChild) .* oriParent;
p2c = project2FundamentalRegion(mean(p2cmori,'robust'));

%% Determine the best fitting plane normal in child coordinates

chi = mean((pi/2-abs(angle_outer(traceChild,n,'min'))).^p,1).^(1/p);

[m,id1] = min(chi(:));

nChild = Miller(n(id1),csChild,'antipodal');
disp(['h_child = ' char(nChild ./ nChild.h)])
disp(['chi = ' xnum2str(m./degree)])

%% Figure 8c) - the traces in child coordinates / corrected

scatter(traceChild)
annotate(vChild,'labeled','backgroundColor','w')

hold on
circle(nChild.symmetrise,'color',ind2color(3),'linewidth',2)
hold off

%saveFigure('../../pic/traceChildCor.png')


%% Figure 8f) - The functional chi^2 in child coordinates / corrected

contourf(Miller(n,csChild),chi ./ degree,'contours',20)
mtexColorMap black2white
mtexColorbar
%setColorRange([min(d(:)),10])

hold on
plot(nChild.symmetrise,'MarkerColor',ind2color(3))
hold off

%saveFigure('../../pic/chiChildCor.png')
