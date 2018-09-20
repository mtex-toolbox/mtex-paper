function ebsd = simEBSD(varargin)

cs = crystalSymmetry;

% spatial coordinates
[opt.y,opt.x] = meshgrid(1:25,1:100);
opt.x = opt.x(:); opt.y = opt.y(:);

f = @(x) pi/10*(1+sign(x-0.4999)).*(1-x);

ori = orientation('Euler',0,f(opt.x./max(opt.x)),0,cs);

ebsd = EBSD(ori,ones(numel(opt.x),1),{cs},'options',opt,'unitCell',....
  calcUnitCell([opt.x,opt.y]));

if check_option(varargin,'Poussin')
  odf = unimodalODF(quaternion.id,'halfwidth',get_option(varargin,'Poussin'));

  ebsd.orientations = discreteSample(odf.components{1},numel(opt.x)) .* ori;
end

if check_option(varargin,'saltPepper')
  odf = uniformODF(cs);

  ind = discretesample(numel(opt.x),numel(opt.x) * get_option(varargin,'saltPepper'));
  ebsd.orientations(ind) = discreteSample(odf,numel(ind));
end

ebsd = ebsd.gridify;

end