% Script file for the publication 
%
%% Approximation of Manifold-valued Functions
% 
% by Ralf Hielscher and Laura Lippert
%
% For running this script MTEX is required
%
%% Section 4.1 Wave Velocities
%

% define the crystal symmetry of Olivine
cs = crystalSymmetry('mmm',[4.7646 10.2296 5.9942],'mineral','Olivin');

% import the single crystal stiffness tensor of Olivine
C = stiffnessTensor.load('Olivine1997PC.GPa',cs);

% include density
rho = 3.355;
C = addOption(C,'density',rho);

% import a distribution function of crystal orientations 
% of a certain polycrystaline specimen
load odf.mat

% compute the polycrystaline stiffness tensor
CV = mean(C,odf)

% 
%T = ChristoffelTensor(C,vector3d.X)


%% Fig. 1a - Polarization Direction of the Fastest Shear Wave

% some spherical grid of propagation directions for plotting
pS = quadratureS2Grid(32);

% compute the polarization direction for these propagation directions 
[~,~,~,~,f] = velocity(CV,pS);

% plot the result
plot(pS,f,'upper','linewidth',2)


%% determine singularities
% the singularities of the vector field of the fastest shear wave appears
% at the positions where the velocities of the fast and slow shear wave
% coincide

[~,vs1,vs2] = velocity(CV,'harmonic','bandwidth',8);

%plot(vs1-vs2)

[~,rSing] = min(vs1-vs2,'numLocal',2)

annotate(rSing,'MarkerSize',15)

%% annotate the input data for the harmonic approximation

% bandwidth of the harmonic approximation
L = 8;

% the Chebyshev quadrature nodes on the sphere for twice the bandwidth
S2 = quadratureS2Grid(2*L);

% compute the polarization direction for these propagation directions 
[~,~,~,~,f] = velocity(CV,S2);

% plot the input data for the Fourier Approximation
hold on
plot(S2,f,'upper','linewidth',3)
hold off

%saveFigure('../pic/ps1Data.pdf','crop')

%% Fig. 1b - Harmonic Approximation 

% perform the harmonic approximation 
[~,~,~,~,f_tilde] = velocity(CV,'harmonic','bandwidth',L);

% plot the result
figure()
plot(f_tilde,'upper','linewidth',2,'color',ind2color(1))

 %saveFigure('../pic/ps1Approx.pdf','crop')

%% Fig. 1c - Error plot

% a very fine grid of propagation directions
r = plotS2Grid('upper');

% the true polarization directions
[~,~,~,~,f,~] = velocity(CV,r);

% the approximated polarization directions
ft_r = f_tilde.eval(r);

% the approximation in the embedding space
f_tilde_rX = f_tilde.sF.eval(r).'; 
f_tilde_rX = f_tilde_rX([1 2 4 2 3 5 4 5 6],:);

% visualize the error |f-I_Mf|
err_func = (sqrt(2-2*dot(ft_r,f.').^2))';

plot(r,err_func,'upper','linewidth',3,'contours',400)

%setColorRange([1e-3,1],'linear')
mtexColorbar('fontSize',26)
mtexColorMap white2black
%mtexColorMap parula

% saveFigure('../pic/ps1error8.pdf','crop')


 %% Fig 1d - The error of the linear approximation

% error without projection
e_f = vecnorm(reshape(double(dyad(f(:),2)),9,[]) - f_tilde_rX);

figure()
plot(r,e_f)
mtexColorbar('fontSize',26)
mtexColorMap white2black

% mark reach by a contour line
hold on
h=plot(r,e_f,'contour',1/sqrt(2),'linewidth',2,'linecolor','b');
hold off

 %saveFigure('../pic/ps1error8.pdf','crop')
 
 
%% Fig. 2a - The derivative Df

% a fine grid of propagation directions
r = plotS2Grid('upper','resolution',0.5*degree);

% determine the derivative numerically,
% r1, r2, are grids on S2 which we use for finite difference

h = 10^(-6);

[t1,t2] = S2VectorField.tangential(r);
r1 = normalize(r + h * t1);
r2 = normalize(r + h * t2);

% point evaluations
[~,~,~,~,f_r] = velocity(CV,r); 
[~,~,~,~,f_d1] = velocity(CV,r1); 
[~,~,~,~,f_d2] = velocity(CV,r2);

% embed into R^3x3 and compute finite difference
t1 = (dyad(f_r(:),2)-dyad(f_d1(:),2))/h;
t2 = (dyad(f_r(:),2)-dyad(f_d2(:),2))/h;

% store as 9x2 matrices
Df = [reshape(double(t1),9,1,[]),reshape(double(t2),9,1,[])];

% norm of the differential of the real function
df_norm = zeros(size(r));
for i = 1:size(Df,3)
  df_norm(i)= norm(Df(:,:,i));
end

figure(1)
plot(r,df_norm,'log','contours',500)
%setColorRange([1e-2,1e4])
%setColorRange([0,1e8])
setColorRange([5e-1,2e2])
cb = mtexColorbar('fontSize',26)
cb.Ticks = [1e0 1e1 1e2 1e3];
mtexColorMap white2black 

% saveFigure('../pic/dps1.png','crop')


%% Fig. 2b - the error of the differential with projection 
%           ||df - d I_M f_tilde||

% we use here a larger polynomial degree
L = 64;

% compute the polarization direction for these propagation directions 
[~,~,~,~,f_tilde] = velocity(CV,'harmonic','bandwidth',L);

f_tilde_r = f_tilde.eval(r);
f_tilde_d1 = f_tilde.eval(r1);
f_tilde_d2 = f_tilde.eval(r2);

% embed into R^3x3 and compute finite difference
t1 = (dyad(f_tilde_r(:),2)-dyad(f_tilde_d1(:),2))/h;
t2 = (dyad(f_tilde_r(:),2)-dyad(f_tilde_d2(:),2))/h;

% store as 9x2 matrices
DI_Mf_tilde = [reshape(double(t1),9,1,[]),reshape(double(t2),9,1,[])];


e = zeros(size(r));
for i = 1:size(DI_Mf_tilde,3)
  e(i)= norm(DI_Mf_tilde(:,:,i)-Df(:,:,i));
end

figure(2)
pcolor(r,e,'upper')
setColorRange([1e-2,1.1])

mtexColorbar('fontSize',26)
mtexColorMap white2black 

%saveFigure('../pic/dps1Error.png','crop')


%% Fig. 2c - The theoretic error bound in Thm 3.2
% 
% This requires some incredients:
%
% (a) the point wise error |f - f_tilde| without projection
% (b) the point wise error |d f - d f_tilde| without projection
% (c) the norm of the derivative |d f| - was computed for Fig. 2a -> dNorm

%% (a) - e_f = |f-f_tilde| without projection

f_tilde_rX = f_tilde.sF.eval(r).'; 
f_tilde_rX = f_tilde_rX([1 2 4 2 3 5 4 5 6],:);

e_f = vecnorm(reshape(double(dyad(f(:),2)),9,[]) - f_tilde_rX);

figure(1)
plot(r,e_f)


%% just for comparison |f-f_tilde| with projection

r = plotS2Grid('upper','resolution',0.5*degree);

[~,~,~,~,f] = velocity(CV,r);

e_fp = norm(dyad(f(:),2) - dyad(f_tilde_r(:),2));

figure(2)
% according to Thm 3.1 thus should be lower than 2*e_f
plot(r,reshape(e_fp,size(r)))

% plot(r,2*e_f(:) - e_fp(:))


%% (b) - edf = |df - df_tilde| without projection

% function values in R^6 for finite differences
f_tilde_rX = f_tilde.sF.eval(r); 
f_tilde_rX1 = f_tilde.sF.eval(r1);
f_tilde_rX2 = f_tilde.sF.eval(r2);

% df_tilde
df_tilde1 = reshape((f_tilde_rX - f_tilde_rX1).' /h, 6,1,[]);
df_tilde2 = reshape((f_tilde_rX - f_tilde_rX2).' /h, 6,1,[]);

% as matrix
Df_tilde = [df_tilde1 df_tilde2];

% expand matrix
Df_tilde = Df_tilde([1 2 4 2 3 5 4 5 6],:,:);

edf_M = Df_tilde - Df;

for i = 1:size(edf_M,3)

  edf(i) = norm(edf_M(:,:,i));

end

%% Fig. 2c - The theoretic error bound in Thm 3.2

% reach 
tau = 1/sqrt(2);    

% right hand side of Theorem 3.2
C = (2/tau + 1./max((tau - e_f),0)) .* (edf + df_norm(:).');
RHS = edf + C .* e_f;


pcolor(r,reshape(RHS,size(r)),'upper')

%setColorRange([1e-2,1e1])
setColorRange([1e-2,1.1])

mtexColorbar('fontSize',26)
mtexColorMap white2black 


%saveFigure('../pic/dps1Bound.png','crop')


%max(e(ind_Omega)./RHS(ind_Omega))
