% MTEX script for the paper
%
%% Habit plane determination from reconstructed parent phase orientation maps
%
% by Frank Niessen, Tuomo Nyyssönen, Azdiar A. Gazder, Ralf Hielscher:
%
% * to run the script at least MTEX 5.9 has to be installed
% 
%% Section 2 - Theorie
% 
%% Data Import 

load('mat/ebsd2.mat')

childGrains = smooth(job.grainsPrior,10);

%% Figure 2(a) traces in specimen coordinates
% The following lines compute traces for all child grains and orders them
% into a table with rows corresponding to the parent grain id and colums
% corresponding to the variants id. Then only those parent grains are
% considered that have childs of all types of variants. Finally, their
% traces are plotted into a spherical projection.

% trace determination
traces = calcTraces(childGrains, [job.mergeId(:), job.variantId(:)]);

% consider only those, for which a parent orientation has been reconstructed
traces = traces(job.isParent,:);

% consider only parent grains that have all 24 variants as childs
has24V = sum(~isnan(traces),2) == 24;
traces = traces(has24V,:);

% plot colorized with respect to parentId
parentId = repmat((1:nnz(has24V)).',1,24);
plot(traces,ind2color(parentId))

%% Figure 2(b) traces in parent coordinates (naive approach)
% In this figure the traces are maped into crystal coordinates using the
% parent grain orientations.

% the corresponding parent orientation
oriParent = job.parentGrains.meanOrientation(has24V);

% turn traces into parent coordinates
tracesParent = inv(oriParent) .* traces;

% plot colorized with respect to parentId
plot(tracesParent,ind2color(parentId))

%% Figure 2(c) traces in parent coordinates including all symmetrically equivalent
% In this figure we also consider all symmetrically equivalent parent
% orientations

plot(tracesParent,ind2color(parentId),'symmetrised')

%% Figure 2(d) traces in parent coordinates using the variant specific parent orientation
% Finally, we consider for each trace the specific parent orientation that
% leads to the variant of the corresponding child grain.


% determine the variant specific parent orientations
oriPVariant = oriParent.project2FundamentalRegion .* ...
  inv(variants(job.p2c)) .* job.p2c;

% turn traces into parent coordinates
tracesParent = inv(oriPVariant) .* traces;

% plot colorized with respect to parentId
plot(tracesParent,ind2color(parentId))
