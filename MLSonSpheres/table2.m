% This code was used to generate table 2 in the paper
% 'An optimal ansatz space for moving least squares approximation on spheres' 

% the spherical fibonacci grids consist of 2*N+1 points
% the goal is to compute the fill distance and separation distance
% and their ratio for fibonacci grids with increasing N 

clear all;

% grid sizes we are going to use for testing 
I = (1:18)';
N = 5 * 2.^I;
sizes = 2 * N + 1;

% initialize vectors for storing parameters of the grids
filldist = 0 * N;
sepdist = 0 * N;

% make a waitbar (for large grids)
w = waitbar(0,'0 %','Name','Computing fill distances');

% compute fill-distances
for l = 1:size(N,1)
    waitbar((l-1)/numel(N), w, sprintf(['Grid %2.0f of %2.0f'], l, numel(N)));
    n = N(l);
    fibgrid = fibonacciS2Grid(2*n+1);
    filldist(l) = fibgrid.filldist;
end
delete(w);

% make a waitbar (for large grids)
w = waitbar(0, '0 %', 'Name', 'Computing separation distances');

% compute separation distances
for l = 1:size(N,1)
    waitbar((l-1)/numel(N), w, sprintf(['Grid %2.0f of %2.0f'], l, numel(N)));
    n = N(l);
    fibgrid = fibonacciS2Grid(2*n+1);
    sepdist(l) = fibgrid.sepdist;
end
delete(w);

% print the results
format long;
[filldist sepdist] / degree
format short;
