classdef detector
  
  properties
    ncols  % width of the detector in pixel
    nrows  % height of the detector in pixel
    dist   % distance detector to scanning surface
    patterCenter  % (x,y) coordinates of the pattern center
    bounds % bounds of the detector in gnonomic projections [xmin, xmax, ymin, ymax]
    x      % coordinates of the pixels in the gnonomic projection
    y      %
  end
  
  properties (Dependent = true)    
    nodesS2 % the pixel positions as points on the sphere
    edges   %
    vertices
  end
  
  properties (Constant)
    proj = gnonomicProjection;
  end
  
  methods
    function det = detector(ncols, nrows, dist, pc)
    
      det.ncols = double(ncols);
      det.nrows = double(nrows);
      det.dist = double(dist);
      det.patterCenter = double(pc);
      
      det.bounds = double([([0,1] - pc(1))*det.ncols/det.nrows,...
        (pc(2) - [1,0])] ./ det.dist) ;
     
      % generate x and y values of the detector positions in the gnonomic projection
      [det.x,det.y] = ndgrid(linspace(det.bounds(1),det.bounds(2),det.ncols),...
        linspace(det.bounds(3),det.bounds(4),det.nrows));
      
    end
    
    function v = get.nodesS2(det)

      % project x,y values onto the sphere
      v = det.proj.iproject(det.x,det.y);
      
    end
        
    
    function e = get.edges(det)
      
      n = det.nodesS2;
      e = [n(1,1),n(1,end),n(end,end),n(end,1)];
            
    end
    
    function v = get.vertices(det)
      
      n = det.nodesS2;
      v = normalize([cross(n(1,1),n(1,end)),cross(n(1,end),n(end,end)),...
        cross(n(end,end),n(end,1)),cross(n(end,1),n(1,1))]);
            
    end
    
    function [x,y] = vec2xy(det,v)
      
      [x,y] = det.proj.project(v);
      
    end
    
    
    function mask = S2CutOffMask(det,delta)
      
      if nargin == 1, delta = 0.1; end
      mask = S2FunHandle(@(v) evalMask(v));
            
      function value = evalMask(v)

        [x,y] = detector.proj.project(v);

        %cutoff = @(t) (1+erf(t ./ delta-2))./2;
        cutoff = @(t) max(0,min(1,t/delta));

        value = cutoff(x - det.bounds(1)) .* cutoff(y-det.bounds(3)) .*  ...
          cutoff(det.bounds(2)-x) .* cutoff(det.bounds(4)-y) .* (v.z > 0);
      end

    end
    
    function pattern = simulatePattern(det,master,ori,flux,bg)
      % simulate a Kikushi pattern given a master pattern
      
      if nargin == 2, ori = rotation.id; end
      pattern = reshape(master.eval(ori \ det.nodesS2),det.nrows,det.ncols);
      
      % add some noise
      if nargin < 4, flux = 0; end
      if nargin < 5, bg = 0; end
      if nargin > 2 && flux>0, pattern = randp(flux*(pattern+bg))-flux*bg; end

      pattern = reshape(pattern,size(det.x)) ;
    end
    
    function [pHarm,plan] = pattern2Fun(det,pattern,varargin)

      plan = getClass(varargin,'struct',struct('mask',[],'S2G',[],'W',[]));
      if isempty(plan.mask)
        plan.mask = det.S2CutOffMask(get_option(varargin,'delta',0.05)); 
      end
      
      if check_option(varargin,'quadrature')
        
        if isempty(plan.S2G)
          [plan.S2G, plan.W] = quadratureS2Grid(2*get_option(varargin,'bandwdith',128));
          plan.maskDiscrete = plan.mask.eval(plan.S2G);
          plan.isInside = plan.maskDiscrete > 0.1;
          plan.maskDiscrete = plan.maskDiscrete(plan.isInside);
        end
        
        % detector positions of the quadrature grid
        [xQ,yQ] = det.vec2xy(plan.S2G(plan.isInside));
        v = zeros(size(plan.S2G));

        v(plan.isInside) = plan.maskDiscrete .* interp2(det.x.',det.y.',pattern.',xQ,yQ);

        pHarm = S2FunHarmonic.quadrature(plan.S2G,v,'weights',plan.W,varargin{:});

      else

        maskDiscrete = plan.mask.eval(det.nodesS2);
        pHarm = S2FunHarmonic.approximation(det.nodesS2(:),maskDiscrete(:) .* pattern(:),'bandwidth',128,varargin{:});
       
      end
    end
  
  end
end