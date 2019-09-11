function peaks = S2PeakDetection(S2F,det,num,delta)

if nargin == 3, delta = -0.12; end

S2R = sphericalRegion(det.edges,delta);
S2F.antipodal = true;

persistent S2G;
if isempty(S2G)
  S2G =equispacedS2Grid('resolution',1.5*degree,'upper');
  S2G(S2R.checkInside(S2G)) = [];
end

[~,peaks] = max(S2F,'numLocal',150,'tolerance',0.5*degree,'startingnodes',S2G,'iterMax',5);
peaks = peaks(~S2R.checkInside(peaks));
[~,~,I] = unique(peaks, 'tolerance', 2*degree);
peaks = normalize(accumarray(I,peaks));

[value,peaks] = max(S2F,'numLocal',num,'tolerance',0.5*degree,'startingnodes',peaks,'iterMax',5,'maxStepSize',1.5*degree);
isInside = ~S2R.checkInside(peaks);
peaks = peaks(isInside);
value = value(isInside);

[~,~,I] = unique(peaks, 'tolerance', 5*degree);

% if peaks are close together take the pake with the highest value
ind = accumarray(I,(1:length(peaks)).',[],@maxValue);

if length(I) > length(ind)
  peaks = peaks(ind);
  [~,peaks] = max(S2F,'numLocal',num,'tolerance',...
    0.5*degree,'startingnodes',peaks,'iterMax',5,'maxStepSize',2*degree);
end

% remove peaks inside the critical region
peaks = peaks(all(angle_outer(det.vertices,peaks,'antipodal')>7*degree));
peaks(S2R.checkInside(peaks)) = [];


function I = maxValue(I)
  
  [~,id] = max(value(I));
  I = I(id);
  
end

end