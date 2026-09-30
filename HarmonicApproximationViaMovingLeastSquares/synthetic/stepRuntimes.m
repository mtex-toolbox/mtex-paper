%% Synthetic Example: Runtime of the Steps of HAMLS
% This script computes Table 1 of the paper. We compare the runtime of
% Step 2 of Algorithm 1, the evaluation of the MLS approximation on the
% Gauss-Legendre quadrature grid, with the runtime of Step 3, the adjoint
% spherical Fourier transform. We consider all bandwidths from 4 to 128,
% $10^4$, $10^5$ and $10^6$ nodes and the polynomial degrees 2, 3 and 4.
%
% All output is written to the folder |results/stepRuntimes|:
%
%  RuntimeMLSQuadrature.txt - Table 1 shows every 16th bandwidth
%  RuntimeMLSQuadrature.mat - individual measurements
%
% The script has to be run from the folder |synthetic| and takes well under
% an hour. The data set of $10^4$ nodes has to coincide with the one written
% by |errorVsBandwidth.m|, which is checked below.

clear
close all

dataDir = 'data';
outDir = fullfile('results','stepRuntimes');
if ~isfolder(outDir), mkdir(outDir); end

bw = 4:128;               % bandwidths
degrees = [2 3 4];        % polynomial degrees of MLS
nNodes = [1e4 1e5 1e6];   % numbers of nodes
nRep = 3;                 % the median over nRep runs is reported

nB = numel(bw);
nN = numel(nNodes);
nD = numel(degrees);

%% The test function, the nodes and the MLS approximations
% The nodes are drawn exactly as in |errorVsBandwidth.m|, with 4 percent of
% them drawn uniformly. The MLS options coincide with those of
% |errorVsBandwidth.m| as well.

uniformFraction = 0.04;
dataSeed = 3;

mlsOpts = {'oF',4,'tangent',true,'regularize',true,'weight','auto', ...
  'use_smooth_delta',true,'use_vor_weights',true};

F = load(fullfile(dataDir,'ToyExampleFun.mat'));  % fun
fun = F.fun * 8 * pi^2 / sum(F.fun);

% the density is the smoothed function, the data are the shifted function
g = smooth(fun,'halfwidth',5*degree);
fun = fun + 300;

nodes = cell(nN,1);
values = cell(nN,1);
for i = 1:nN
  m = round(uniformFraction * nNodes(i));
  rng(dataSeed)
  nodes{i} = [discreteSample(g,nNodes(i) - m); vector3d.rand(m)];
  values{i} = fun.eval(nodes{i});
end

% check the 1e4 data set against the one of errorVsBandwidth.m
refFile = fullfile('results','errorVsBandwidth','ToyExampleData.mat');
refData = load(refFile);  % fun, nodes, values
dNodes = max(angle(nodes{1},refData.nodes) / degree);
dVals = max(abs(values{1} - refData.values));
if dNodes > 1e-9 || dVals > 1e-9
  warning('The data set differs from %s (%.3g degrees, %.3g in value).', ...
    refFile,dNodes,dVals);
end

mls = cell(nN,nD);
for i = 1:nN
  for j = 1:nD
    fprintf('MLS construction: %.0e nodes, degree %i\n',nNodes(i),degrees(j));
    mls{i,j} = S2FunMLS(nodes{i},values{i},'degree',degrees(j),mlsOpts{:});
  end
end

%% Runtime of Step 2 and Step 3
% The repetitions form the outer loop, so that a drift of the machine over
% time spreads over all bandwidths. Step 3 is timed once per bandwidth and
% repetition, with the values of the first data set and degree, since its
% cost depends only on the grid.

tMLS = nan(nB,nN,nD,nRep);  % Step 2
tQuad = nan(nB,nRep);       % Step 3

% warm up, as the first call carries JIT and allocation overhead
S2G = quadratureS2Grid(bw(1),'GaussLegendre');
for i = 1:nN
  for j = 1:nD
    v = mls{i,j}.eval(S2G);
  end
end
S2FunHarmonic.adjoint(S2G,v);

for r = 1:nRep
  for n = 1:nB
    fprintf('repetition %i of %i, bandwidth %i\n',r,nRep,bw(n));
    S2G = quadratureS2Grid(bw(n),'GaussLegendre');

    for i = 1:nN
      for j = 1:nD
        tic
        v = mls{i,j}.eval(S2G);
        tMLS(n,i,j,r) = toc;
        if i == 1 && j == 1, vQuad = v; end
      end
    end

    tic
    S2FunHarmonic.adjoint(S2G,vQuad);
    tQuad(n,r) = toc;
  end
end

% reshape instead of squeeze keeps the order (bandwidth, nodes, degree)
tMLSmed = reshape(median(tMLS,4,'omitnan'),nB,nN,nD);
tQuadmed = median(tQuad,2,'omitnan');

save(fullfile(outDir,'RuntimeMLSQuadrature'),'bw','nNodes','degrees', ...
  'tMLS','tQuad','tMLSmed','tQuadmed');

%% Write the table
% Table 1 is typeset from this file with pgfplotstable.

writeTable(fullfile(outDir,'RuntimeMLSQuadrature.txt'), ...
  [{'bandwidth','Quadrature'},mlsColumnNames(nNodes,degrees)], ...
  [bw', tQuadmed, reshape(tMLSmed,nB,nN*nD)]);

%% Summary
% We print the numbers quoted in Section 3.2, in particular the ratio of the
% runtimes of Step 2 for degree 4 and degree 2 at the rows of Table 1.

nLast = nB;  % L = 128
jK2 = find(degrees == 2,1);
jK4 = find(degrees == 4,1);
MQ = numel(quadratureS2Grid(bw(nLast),'GaussLegendre'));
step2 = reshape(tMLSmed(nLast,:,jK4),1,[]);

nTable = find(mod(bw,16) == 0);
degFactor = tMLSmed(nTable,:,jK4) ./ tMLSmed(nTable,:,jK2);

fprintf('\nSection 3.2, Table 1\n');
fprintf('L = %i: %i quadrature nodes\n',bw(nLast),MQ);
fprintf('L = %i: Step 3 %.3f s, Step 2 (degree 4) %.2f to %.2f s\n', ...
  bw(nLast),tQuadmed(nLast),min(step2),max(step2));
fprintf('L = %i: Step 2 is %.0f times slower than Step 3 for N = %.0e\n', ...
  bw(nLast),tMLSmed(nLast,1,jK4) / tQuadmed(nLast),nNodes(1));
fprintf('\nStep 2, degree 4 over degree 2\n');
fprintf('    L %s\n',sprintf('   N = %.0e',nNodes));
for t = 1:numel(nTable)
  fprintf('  %3i %s\n',bw(nTable(t)),sprintf('%12.2f',degFactor(t,:)));
end

paperNumbers = struct('MQ128',MQ,'tQuad128',tQuadmed(nLast), ...
  'step2Range128',[min(step2), max(step2)],'degFactor',degFactor, ...
  'stepFactor128',tMLSmed(nLast,1,jK4) / tQuadmed(nLast));

%% A quick look at the results

figure(4); clf
semilogy(bw,tQuadmed,'b-','LineWidth',1.5); hold on
styles = {':','-.','--'};
colors = {[0 0.6 0],[1 0.5 0],[1 0 0]};
for j = 1:nD
  for i = 1:nN
    semilogy(bw,tMLSmed(:,i,j),styles{i},'Color',colors{j});
  end
end
hold off
xlabel('bandwidth L'); ylabel('runtime in seconds'); grid on
title('Step 2 for degree 2, 3, 4 (green, orange, red) and Step 3 (blue)');

%% Helper functions

function names = mlsColumnNames(nNodes,degrees)
% column names MLS_1e<k>_deg<d>, the number of nodes varying fastest

names = cell(1,numel(nNodes)*numel(degrees));
c = 0;
for j = 1:numel(degrees)
  for i = 1:numel(nNodes)
    c = c + 1;
    names{c} = sprintf('MLS_1e%i_deg%i',round(log10(nNodes(i))),degrees(j));
  end
end

end


function writeTable(fname,header,data)
% write a tab separated table with a single header line

assert(numel(header) == size(data,2), ...
  '%s: %i column names for %i columns',fname,numel(header),size(data,2));

fid = fopen(fname,'w');
fprintf(fid,'%s',header{1});
fprintf(fid,'\t%s',header{2:end});
fprintf(fid,'\n');
fclose(fid);

writematrix(data,fname,'Delimiter','tab','WriteMode','append');

end
