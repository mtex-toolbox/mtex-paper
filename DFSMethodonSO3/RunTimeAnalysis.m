%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        RunTimeAnalysis of the NFSOFT-Algorithms (WCT,NFFT,FFT) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFTW 3.3.10, NFFT 3.5.3, and MTEX 6.1


% LEGENDE: 
%   (a) fast WCT w.r.t Chebychev-Polynomials
%   (b) direct WCT
%   (c) NFFT

clear all
% number of runs to average the time
maxIter = 100;

% bandwidth
bandwidth = [2.^(0:8),448,2^9];
ind = 0:length(bandwidth)-1;

% results
Time = struct('bandwidth',bandwidth,'time',NaN(length(ind),13,maxIter), ...
  'columns',["PrecomputationsFPT","fastWCT_withoutPre","fastWCT",...
             "directWCT","NFFT3_N3_m4_sig15"]);

for iter=1:maxIter

  % Save some time (do not iterate to often for high bandwidths)
  if iter>20
    bandwidth(bandwidth>128) = [];
    ind = 0:length(bandwidth)-1;
  end

  for n=ind

    N = bandwidth(n+1);

    % nfft parameters
    m = 4;
    sigma = 1.5;

    % harmonic coefficients
    fhat = rand(deg2dim(N+1),2)*[2;2i] -1-1i;
    
    % nodes for nfft
    M = deg2dim(N+1);
    rot = Euler(rotation.rand(M),'nfft').';


    % (a) fast WCT
    if  (N<=128) || (N>128 && N<=256 && iter<=10)
      tic
      % flags
      % 2^0  -> L_2 normalized Wigner-D functions
      % 2^2  -> use internally usually slower direct DPT algorithm in favor of the fast FPT algorithm 
      % 2^4  -> Wigner-D functions will be normalized such that they satisfy the representation property of the spherical harmonics as defined in the NFFT software package
      % 2^20 -> Do not use the nfft. Only do the Wigner transform and transform the harmonic coefficients to Fourier coefficients.
      nfsoft_flags = bitor(2^4,2^0)+2^20;
      % initialize nfsoft plan
      plan = nfsoftmex('init',N,1,nfsoft_flags,0,1,1000,2*N+2);
      s = toc;
      % set Fourier coefficients
      nfsoftmex('set_f_hat',plan,fhat);
      % fast SO(3) fourier transform
      nfsoftmex('trafo',plan);
      % get function values from plan and normalize
      ghata = nfsoftmex('get_g_hat',plan)*sqrt(8*pi^2);
      t = toc;

      clear ghata

      Time.time(n+1,1,iter) = s;
      Time.time(n+1,2,iter) = t-s;
      Time.time(n+1,3,iter) = t;
      

      % kill plan
      nfsoftmex('finalize',plan);
    end


    % (b) direct WCT
      % flags: 2^0 -> use L_2-normalized Wigner-D functions
      %        2^1 -> make size of result even
      %        2^4 -> use right and left symmetry
  
      F = SO3FunHarmonic(fhat);
      tic
      flags = 2^0+2^1;
      ghatc = wignerTrafo(F,flags,'bandwidth',N);
      t = toc;

      Time.time(n+1,4,iter) = t;


    % (c) NFFT
      ghat = ghatc(:);
      clear ghatc
    if (N<=256) || (N>256 && N<=450 && iter<=10)
      tic
      % initialize nfft plan
      NN = 2*N+2;
      FN = 2*ceil(sigma/2*NN);
      % {FFTW_ESTIMATE} or 64 - Specifies that, instead of actual measurements of different algorithms,
      %                         a simple heuristic is used to pick a (probably sub-optimal) plan quickly.
      %                         It is the default value
      % {FFTW_MEASURE} or 0   - tells FFTW to find an optimized plan by actually computing several FFTs and
      %                         measuring their execution time. This can take some time (often a few seconds).
      fftw_flag = int8(64);
      plan = nfftmex('init_guru',{3,NN,NN,NN,M,FN,FN,FN,3,int8(0),fftw_flag});
      % set rotations as nodes in plan
      nfftmex('set_x',plan,rot);
      % node-dependent precomputation
      nfftmex('precompute_psi',plan);
      % set Fourier coefficients
      nfftmex('set_f_hat',plan,ghat);
      % fast fourier transform
      nfftmex('trafo',plan);
      % get function values from plan
      fc = nfftmex('get_f',plan);
      t = toc;

      Time.time(n+1,5,iter) = t;

      % kill plan
      nfftmex('finalize',plan);
      clear fc
    end


  end      
  progress(iter,maxIter)
end





%% Compute Mean Values

sz = size(Time.time);
A = zeros(sz(1),sz(2));
for k=1:sz(1)
  for l=1:sz(2)
    v = sort(squeeze(Time.time(k,l,:)));
    v = v(~isnan(v));
    if isempty(v)
      A(k,l) = NaN;
    else
      A(k,l) = mean(v(ceil(0.25*length(v)):ceil(0.75*length(v))),'omitnan');
    end
  end
end
% Only use bw = 448 if bw = 512 is not computable 
A(end-1,~isnan(A(end,:))) = NaN;



%% Plot Data

bw = Time.bandwidth;

loglog(bw,A(:,4), 'DisplayName', 'Direct Wigner Transform');
hold on
loglog(bw,A(:,3), 'DisplayName', 'Wigner Transform via FPT with Precomputations');
loglog(bw,A(:,2), 'DisplayName', 'Wigner Transform via FPT without Precomputations');
loglog(bw,A(:,5), 'DisplayName', 'NFFT');
hold off

grid on;
box on;

% x-axis: powers of 2
xticks(2.^(0:9));   % 1,2,4,...,512
xlim([1 512]);

% y-axis: powers of 10
yticks(10.^(-5:2)); % 1e-5,...,1e2
ylim([8e-6 500]);

xlabel('Bandwidth N');
ylabel('CPU time (s)');

legend('Location', 'northwest');

