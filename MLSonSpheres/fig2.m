
% Moving Least Squares for different bases on Fibonacci-Grids
clear all;

% precomputed fill distances in degree
filldist = [
  33.904174188722251
  24.315997333536416
  17.276216676397190
  12.160594123149252
  8.564595202533692
  5.981372838649667
  4.266094352447787
  2.885579634060381
  2.150696540997295
  1.432226735679432
  1.002601683475767
  0.700404799231414
  0.494254792603826 
  0.352863964441357
  0.248635554143566
  0.176626077525876
  0.123640549560582
  0.087185772736235
  0.061867254619169
  0.043650384248674
  0.031067484657432
  0.021938527395164
  0.015655710838108
  0.011038833120364
  0.007674297907596
  0.005442403665982];

% polynomial degree and dimension of the ansatz spaces
deg = 3;
m_all = (deg+1)^2;
m_even = (deg+1) * (deg+2) / 2;

% generate random test points
rng(0);
s = 1000;
v = vector3d.rand(s,1);

% evaluate f on T for error analysis
tv = testfun(v);

% initialize arrays for the loop over all grids
I = (1:18)';
M = numel(I);
N = 5 * 2.^I;
filldist = filldist(I);
sizes = 2*N + 1;

% arrays for storing data
% pointwise errors for each fill distance
err_all_harm = zeros(s,numel(I));
err_even_harm = zeros(s,numel(I));
err_even_mon_cent = zeros(s,numel(I));
err_tangent = zeros(s,numel(I));

% pointwise conditions for each fill distance
cond_all_harm = zeros(s, numel(I));
cond_even_harm = zeros(s, numel(I));
cond_even_mon_cent = zeros(s, numel(I));
cond_tangent = zeros(s, numel(I));

% choose the support radius delta
deltascale_all  = 4.5;
deltascale_even = 3.5;
delta_all  = deltascale_all  * filldist * degree;
delta_even = deltascale_even * filldist * degree;

% create a waitbar (for running on large grids)
w = waitbar(0, '0 %', 'Name', 'MLS on Fibonacci-Grids');

% loop over all grid sizes and compute the MLS approximation on the test nodes
for i = 1 : M
  % progress in the loop
  p = (i - 1) / M;

  % make the grid and evaluate f on it
  n = N(i);
  fg = fibonacciS2Grid(2 * n + 1);
  f = testfun(fg);

  % all harmonics 
  waitbar(p, w, sprintf('All harmonics (%2.0f of %2.0f)', i, M));
  sF = S2FunMLS(fg, f, deg, delta_all(i), 'all_degrees', ...
    'weight', 'squared hat');
  [vals, cond_all_harm(:,i)] = sF.eval(v);
  err_all_harm(:,i) = abs(tv - vals);
  
  % even harmonics
  waitbar(p, w, sprintf('Even harmonics (%2.0f of %2.0f)', i, M));
  sF = S2FunMLS(fg, f, deg, delta_even(i), 'weight', 'squared hat');
  [vals, cond_even_harm(:,i)] = sF.eval(v);
  err_even_harm(:,i) = abs(tv - vals);

  % even monomials centered
  waitbar(p, w, sprintf('Even monomials centered (%2.0f of %2.0f)', i, M));
  sF = S2FunMLS(fg, f, deg, delta_even(i), 'monomials', 'centered',... 
      'weight', 'squared hat');
  [vals, cond_even_mon_cent(:,i)] = sF.eval(v);
  err_even_mon_cent(:,i) = abs(tv - vals);

  % pullbacks of polynomials to the tangent space
  waitbar(p, w, sprintf('Tangent Space (%2.0f of %2.0f)', i, M));
  sF = S2FunMLS(fg, f, deg, delta_even(i), 'tangent',... 
      'weight', 'squared hat');
  [vals, cond_tangent(:,i)] = sF.eval(v);
  err_tangent(:,i) = abs(tv - vals);
end

% maximal errors among all centers, one value for every fill distance
maxerr_all_harm = max(err_all_harm, [], 1);
maxerr_even_harm = max(err_even_harm, [], 1);
maxerr_even_mon_cent = max(err_even_mon_cent, [], 1);
maxerr_tangent = max(err_tangent, [], 1);

% compute the maximal condition number of the gram matrix
% for each fill distance for each method
maxcond_all_harm = max(cond_all_harm, [], 1);
maxcond_even_harm = max(cond_even_harm, [], 1);
maxcond_even_mon_cent = max(cond_even_mon_cent, [], 1);
maxcond_tangent = max(cond_tangent, [], 1);

% assemble the max conds for each method into one array
maxconds = [maxcond_all_harm; maxcond_even_harm; ...
  maxcond_even_mon_cent; maxcond_tangent]';
% same for errors
errors = [maxerr_all_harm; maxerr_even_harm; ...
  maxerr_even_mon_cent; maxerr_tangent];

% delete the wait bar
delete(w);

% plot the results
% 1 - errors
figure(1);
loglog(filldist, errors, 'LineWidth', 2);
title('Error vs. fill distance');
[~, hobj, ~, ~] = ...
  legend({'All harmonics', 'Even harmonics', ... 
  'Even monomials centered', 'Tangent Space'}, ...
  'Location','southeast', 'FontSize', 20);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',4);

% 2 - conditions of the gram matrices
figure(2);
loglog(filldist, maxconds, 'LineWidth', 2);
title('Condition vs. fill distance');
[~, hobj, ~, ~] = ...
  legend({'All harmonics', 'Even harmonics', ... 
  'Even monomials centered', 'Tangent Space'}, ...
  'Location','southeast', 'FontSize', 20);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',4);

% the test function
function vals = testfun(v)
% s = size(v, 1);
n = [1 1 2 1 1]';   
alpha = [5 7 6 5 2.1]';   
c = [2 0.5 -2 -2 0.2];
centers = [       0         0  1       ;
           0.932039         0  0.362358;
          -0.362154  0.612280  0.696707;
           0.904035  0.279651 -0.323290;
         -0.0479317 -0.424684 -0.904072];
% vals = (c * exp(-repmat(alpha,1,s) .* (1 - centers * v.xyz').^repmat(n,1,s)))';
vals = exp(v.x .* v.y .* v.z);
end