% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      SIMULATION OF THE ERROR OF THE WIGNER COEFFICIENT TRANSFORM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFTW 3.3.10, NFFT 3.5.3, and MTEX 6.1


%   (a) direct WCT
%   (b) fast WCT by usage of the FPT with stabilization step
%   (c) fast NFSOFT (fast WCT [with stabilization step] and NFFT)
%   (d) fast WCT by usage of the FPT without stabilization step

clear all

% nfsoftmex-Flags
% 2^0  -> L_2 normalized Wigner-D functions
% 2^2  -> use internaly usually slower direct DPT algorithm in favor of the fast FPT algorithm
% 2^4  -> Wigner-D functions will be normed such that they satisfy the representation property of the spherical harmonics as defined in the NFFT software package
% 2^13 -> If this flag is set, the fast NFSOFT algorithms (see nfsoft_trafo, nfsoft_adjoint) will use internally the FPT algorithm without the stabilization scheme and thus making bigger errors for higher bandwidth but becoming significantly faster
% 2^20 -> Do not use the nfft. Only do the Wigner transform and transform the harmonic coefficients to Fourier coefficients.

m = 4;
sigma = 3;
kappa = 10^3;

bw = [1,2,4,8,16,20,25,32,35:5:60,64,70:5:120,128,130:10:200,225,256];
maxiter = 100;
for iter = 1:maxiter
  for N = bw
    fprintf([num2str(N),' '])
    for s = [0,1/2,1,3/2]
      
      % harmonic coefficients
        % fhat = rand(deg2dim(N+1),2)*[2;2i]-1-1i;
        fhat = sqrt(rand(deg2dim(N+1),1)) .* exp(2i*pi*rand(deg2dim(N+1),1));
        if s ~= 0
          for n=0:N
            ind = deg2dim(n)+1:deg2dim(n+1);
            fhat(ind) = fhat(ind) / norm(fhat(ind)) / (2*n+1)^(s);
          end
          fhat = fhat/sqrt(N+1);        
        end

      % (a) Direct WCT

        ghat = wignerTrafo(SO3FunHarmonic(fhat),2^0+2^4,'bandwidth',N);
        ghat = Weights(ghat,N);
        ghat = (1i).^(reshape(-N:N,1,1,[]) - (-N:N).') .* ghat;
        fhatWCT = wignerTrafoAdjointmex(N,ghat,2^0,[1,1,1,1]);

      % (b) Fast Polynom Transform with FFT

        % forward FPT
        nfsoft_flags = bitor(2^4,2^0)+2^20;
        plan = nfsoftmex('init',N,1,nfsoft_flags,0,m,kappa,2*ceil(sigma/2*(2*N+2)));
        nfsoftmex('set_f_hat',plan,fhat);
        nfsoftmex('trafo',plan);
        ghat = nfsoftmex('get_g_hat',plan)*sqrt(8*pi^2);

        % Weighting
        ghat = reshape(ghat,2*N+2,2*N+2,2*N+2);
        ghat = permute(ghat,[3,1,2]);
        ghat(2:end,2:end,2:end) = Weights(ghat(2:end,2:end,2:end),N);
        ghat = permute(ghat,[2,3,1]);

        % adjoint DPT
        nfsoftmex('set_g_hat', plan, ghat(:));
        nfsoftmex('adjoint', plan);
        fhatFPT = nfsoftmex('get_f_hat', plan)*(sqrt(8)*pi);
        nfsoftmex('finalize',plan);

      % (c) FPT with NFSOFT instead of FFT
        nfsoft_flags = bitor(2^4,1);
        r = quadratureSO3Grid(N);
        plan = nfsoftmex('init',N,length(r),nfsoft_flags,0,m,kappa,2*ceil(sigma/2*(2*N+2)));
        nfsoftmex('set_x',plan,Euler(r,'nfft').');
        nfsoftmex('precompute',plan);
        nfsoftmex('set_f_hat',plan,fhat);
        nfsoftmex('trafo',plan);
        v = nfsoftmex('get_f',plan) * (sqrt(8)*pi);
        nfsoftmex('finalize',plan);
        fhatNFSOFT = calcFourier(SO3FunHarmonic.quadrature(r,v,'nfsoft'));
        clear r      
        
      % (d) Fast Polynom Transform without Stabilization step

        % forward FPT
        nfsoft_flags = bitor(2^4,2^0)+2^20+2^13;
        plan = nfsoftmex('init',N,1,nfsoft_flags,0,m,kappa,2*ceil(sigma/2*(2*N+2)));
        nfsoftmex('set_f_hat',plan,fhat);
        nfsoftmex('trafo',plan);
        ghat = nfsoftmex('get_g_hat',plan)*sqrt(8*pi^2);

        % Weighting
        ghat = reshape(ghat,2*N+2,2*N+2,2*N+2);
        ghat = permute(ghat,[3,1,2]);
        ghat(2:end,2:end,2:end) = Weights(ghat(2:end,2:end,2:end),N);
        ghat = permute(ghat,[2,3,1]);

        % adjoint DPT
        nfsoftmex('set_g_hat', plan, ghat(:));
        nfsoftmex('adjoint', plan);
        fhatFPTwos = nfsoftmex('get_f_hat', plan)*(sqrt(8)*pi);
        nfsoftmex('finalize',plan);

      % Error estimation

        load('ErrorWCT.mat')
        e_l2 = [norm(fhat-fhatWCT),norm(fhat-fhatFPT),norm(fhat-fhatNFSOFT),norm(fhat-fhatFPTwos)];
        E(end+1,:) = [N,e_l2/norm(fhat,1)];
        save('ErrorWCT','E','columns')
    end
  end
  progress(iter,maxiter)
end


%% Function

function ghatH = Weights(ghat,N)
  
  % Evaluation
  val = fftn(ghat,[2*N+2,4*N,2*N+2]);
  val = val(:,1:2*N+1,:);
  z = (0:2*N+1).' * (N/(N+1)) + 0.5*(0:2*N) + reshape(0:2*N+1,1,1,[]) * (N/(N+1));
  val = exp(1i*pi*z) .* val;
  
  % Quadrature
  W = quadratureSO3Grid.fclencurt2(2*N+1)'/(8*(N+1)^2);
  ghatH = ifftn( W.* val ,[2*N+2,4*N,2*N+2]);
  ghatH = ifftshift(ghatH);
  ghatH = 16*N*(N+1)^2 * ghatH(2:end,N+1:3*N+1,2:end,:);

end






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        ERROR-ANALYSIS OF THE WIGNER COEFFICIENT TRANSFORM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Preprocessing

% Compute Mean and Variance of E over same bandwidth
bw = E(:,1);
E = E(:,2:end);
A = []; V=[];
for n = unique(bw')
  e = E(bw==n,:);
  A(end+1,:) = mean(e,1);
  V(end+1,:) = var(e,0,1);
end
bw = unique(bw);

% For the FWT-time, we use either (b) or (c), depending on which is faster.
% [directWCT, FWT, FWTwoStab]
M = [A(:,1),min(A(:,[2,3]),[],2),A(:,4)];
V = [V(:,1),min(V(:,[2,3]),[],2),V(:,4)];



%% Plot the Error

loglog(bw,M(:,1),'DisplayName','Direct Wigner Transform')
hold on
loglog(bw,M(:,2),'DisplayName','Wigner Transform via FPT')
loglog(bw,M(:,3),'DisplayName','Direct Wigner Transform via FPT without stab. step')
hold off

grid on;
box on;

% x-axis: powers of 2
xticks(2.^(0:8));   % 1,2,4,...,256
xlim([1 256]);

% y-axis: powers of 10
yticks(10.^(-16:4:0));
ylim([1e-18 100]);

xlabel('Bandwidth N');
ylabel('$E_{\ell_1 \to \ell_2}$ Error', 'Interpreter', 'latex');

legend('Location', 'northwest');




%% Plot the Variance

loglog(bw,V(:,1),'DisplayName','Direct Wigner Transform')
hold on
loglog(bw,V(:,2),'DisplayName','Wigner Transform via FPT')
loglog(bw,V(:,3),'DisplayName','Direct Wigner Transform via FPT without stab. step')
hold off

grid on;
box on;

% x-axis: powers of 2
xticks(2.^(0:8));   % 1,2,4,...,256
xlim([1 256]);

% y-axis: powers of 10
yticks(10.^(-45:10:0));
ylim([1e-45 1]);

xlabel('Bandwidth N');
ylabel('Variance of the $E_{\ell_1 \to \ell_2}$ Error', 'Interpreter', 'latex');

legend('Location', 'northwest');

