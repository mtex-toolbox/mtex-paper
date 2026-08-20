function [Erel,Eabs] = discrepancySO3(f,ori,c,varargin)
% kernel discrepancy between a weighted orientation set and an ODF
%
% Description
%
% This is the functional that <SO3Fun.optimalSample.html |optimalSample|>
% minimizes, i.e. the squared distance of the discrete measure
%
%   mu = lambda * sum_j c_j delta_{R_j},   lambda = int_{SO(3)} f,
%
% and the ODF f in the norm of the restricted distance kernel,
%
%   J(R,c) = sum_{l=1}^{N} 8*pi^2*A_l/(2l+1) * sum_{k,k'} |muhat - fhat|^2 .
%
% The degree 0 term is dropped, since the kernel is only conditionally
% positive definite - it vanishes anyway as long as the weights sum up to 1,
% see |optimalSample|. In contrast to the L2 error of a kernel density
% estimate, this measure does not depend on any smoothing parameter and it is
% the quantity the sampling algorithm actually optimizes: it bounds the
% quadrature error of the orientation set for every property of bandwidth N.
%
% The absolute value J is normalized by the same norm of f itself, so that
% Erel = sqrt(J)/||f||_psi is dimensionless and comparable across sample
% sizes and textures.
%
% Syntax
%   Erel = discrepancySO3(f,ori)             % equal weights 1/M
%   Erel = discrepancySO3(f,ori,c)
%   [Erel,Eabs] = discrepancySO3(f,ori,c,'bandwidth',32)
%
% Input
%  f   - @SO3Fun, ODF
%  ori - @orientation, @rotation
%  c   - weights (non negative, are normalized to sum 1; default 1/M)
%
% Output
%  Erel - relative discrepancy sqrt(J)/||f||_psi
%  Eabs - sqrt(J)
%
% Options
%  bandwidth - harmonic degree taken into account (default = 32), has to be
%              the bandwidth the sample was optimized for
%
% See also
% SO3Fun/optimalSample SO3RestrictedDistanceKernel discrepancyS2

bw = get_option(varargin,'bandwidth',32);

f = SO3FunHarmonic(f,'bandwidth',bw);
f.bandwidth = bw;

M = numel(ori);

if nargin < 3 || isempty(c), c = ones(M,1)/M; end
c = c(:);
c = c/sum(c);

% the nodes have to carry the symmetries of f, since the difference of the
% discrete measure and f is formed below
ori = orientation(ori(:),f.CS,f.SS);

% Chebyshev coefficients of the kernel and the resulting norm weights, zero
% in degree 0
psi = SO3RestrictedDistanceKernel(bw+1);

w = zeros(deg2dim(bw+1),1);
for l = 1:bw
  w(deg2dim(l)+1:deg2dim(l+1)) = sqrt( 8*pi^2 * psi.A(l+1)/(2*l+1) );
end

% Wigner coefficients of the discrete measure mu. Note the scaling by
% sqrt(8)*pi, which converts between the normalization of MTEX and the
% orthonormal Wigner-D functions, exactly as in optimalSample.
lambda = sum(f);
mu = SO3FunHarmonic.adjointNFSOFT(ori,c,'bandwidth',bw);
mu.bandwidth = bw;

D = (lambda/(sqrt(8)*pi)) * mu - (sqrt(8)*pi) * f;
D.bandwidth = bw;

Eabs = sqrt(sum(abs(w.*D.fhat).^2));
Erel = Eabs ./ sqrt(sum(abs(w.*((sqrt(8)*pi)*f.fhat)).^2));

end
