%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Anisotropic Rotational Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run with MatlabR2025b


clear
close all

D = diag([0.15,0.4,1]);
beta = 2;
eta = 0.3;

n0 = zvector;
Q = rotation.byAxisAngle(yvector,70*degree);
n1 = Q * n0;

%% Potential

function u = Ufun(r,n)
  s = size(r);
  r = r(:);
  eta = 0.3;
  cz = dot(n,r*zvector);
  cx = dot(n,r*xvector);
  cy = dot(n,r*yvector);
  u = -0.5 .* (3 .* cz.^2 - 1) -eta .* (cx.^2 - cy.^2);
  u = reshape(u,s);
end

%% Initial equilibrium distribution

f0 = SO3FunHandle(@(r) exp(-beta .* Ufun(r,n0)));
f0 = SO3FunHarmonic(f0);

% normalize as ODF
f0 = f0 ./ sum(f0);

%% Potential after director switching

U = SO3FunHandle(@(r) Ufun(r,n1));
U = SO3FunHarmonic(U);


%% Simulation


numIter = 3000;
L = 32;

% precompute
gradU = grad(U,'right');

f = f0;

for k=1:numIter
  fprintf([num2str(k),' / ',num2str(numIter),'  '])

  G = grad(f,'right') + (beta*f).*gradU;
  G.SO3F = D*G.SO3F;
  f = f + 1/numIter * div(G);
  if f.bandwidth>L, f.bandwidth = L; end

  % Plot
  if mod(k,10)==0
    figure(1)
    plot3d(f,'AxisAngle')
    setColorRange([0.0016,0.0584])
    h = gcf();
    set(h,'position',[10,10,1000,1000])
    drawnow
    exportgraphics(h,['PlotAxisAngle3d/AxisAngle3dPic',num2str(k),'.png'],'Resolution',300,'Padding',100);
    drawnow
    close(h)
  end

end

%%

figure(1)
plot3d(f0,'AxisAngle')
setColorRange([0.0016,0.0584])
h = gcf();
set(h,'position',[10,10,1000,1000])
drawnow
exportgraphics(h,'PlotAxisAngle3d/AxisAngle3dPic0.png','Resolution',300,'Padding',100);
drawnow
close(h)




%%

finf = SO3FunHandle(@(r) exp(-beta .* Ufun(r,n1)));
finf = finf / sum(finf) * 8*pi^2

figure(3)
plot3d(finf,'AxisAngle')

%%

