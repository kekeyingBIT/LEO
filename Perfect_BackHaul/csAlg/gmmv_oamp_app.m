function [X_est_mu, lambda,N_iteration] = gmmv_oamp_app(R,Phi,damp,niter,nvar,X,Supp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% Y: Measurements
% S_wave: Meansurement matrix
% damp: Damping value
% Iiter: The maximum number of iterations
% prior: 0/1 for Gaussian/MPSK prior
% L: Modulation order of MPSK
% Output:
% miu_storage：The posterior mean of each iteration
% lambda: The posterior sparsity ratio
% thegma2; The estimated noise variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse Dimension 
[M,T,P] = size(R);
[~,N,~] = size(Phi);

    
%% initialize for EM parameter

% initialize noise
% sigma2 = norm(R(:),2)^2/((1+10^(SNRdB/10)*M*T*P)); %
sigma2 = nvar;
% initialize sparsity ratio
alpha = M/N;
normal_cdf = @(x) 1/2.*(1+erf(x/sqrt(2)));
normal_pdf = @(x) 1/sqrt(2*pi).*exp(-x.^2/2);
alpha_grid = linspace(0,10,1024);
rho_SE = (1 - 2/alpha.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)))...
    ./ (1 + alpha_grid.^2 - 2.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)));
rho = alpha*max(rho_SE);
% initialize parameter mean
mux = 0;
% initialize parameter variance, (realiable when zero mean, ||A||_F^2*rho*vx2*T + M*T*sigma2 = || Y-A*x_hat0||_F^2)
vx2 = ( ( norm(R(:),2)^2-T*M*P*sigma2 ) / (T*norm(Phi(:),2)^2*rho) );

% sigma2 = sigma_true;
% rho = rho_true;
% mux = mux_true;
% vx2 = vx2_true;

%% initialize oamp algorithm

X_hat = zeros(N,T,P);
V2 = vx2*rho*ones(P,1);

X_hat_pre = X_hat;
V2_pre = V2;

Sigma2 = sigma2*ones(P,1);
Mux = mux*ones(P,1);
Vx2 = vx2*ones(P,1);
Rho = rho*ones(N,T,P);

X_est_mu = zeros(N,T,P,niter); % store result for output
X_hat_post = zeros(N,T,P); % store result for cal NMSE
MSE_simulation = zeros(1,niter);  % store result for visualization

% precompute W_t for acceralate
% W_t_P = zeros(N,M,P);
% TrBt = zeros(P,1);
% TrWt = zeros(P,1);
% for p = 1:P
%     A = Phi(:,:,p);
%     W_hat = A'/(A*A');
%     W_t = N/trace(W_hat*A)*W_hat;
%     B_t = eye(N)-W_t*A;
%     W_t_P(:,:,p) = W_t;
%     TrBt(p) = 1/N*trace(B_t*B_t');
%     TrWt(p) = 1/N*trace(W_t*W_t');
% end

for it = 1:niter
    for p = 1:P  % estimate different pages separately
        %% parse value and parameter
        A = Phi(:,:,p);
        Y = R(:,:,p);

        H_hat = X_hat(:,:,p);
        v2 = V2(p);

        H_hat_pre = X_hat_pre(:,:,p);
        v2_pre = V2_pre(p);

        sigma2 = Sigma2(p);
        rho = Rho(:,:,p);
        mux = Mux(p);
        vx2 = Vx2(p);
        %% LE: LMMSE estimator

%         method 1  (39a)
%         W_hat = v2*A' / ( v2*(A*A')+sigma2*eye(M) );
%         W_t = K/trace(W_hat*A)*W_hat;  % Decorrelated
%         x_post = x_hat + W_t*(y-A*x_hat);  % posterior mean
%         V_post = v2*( eye(K) - v2*A'* ( v2 * (A*A') + sigma2 * eye(M) )^(-1) * A);  % posterior variance
%         v_post = trace(real(V_post))/K;
%         v_ext = 1/(1/v_post - 1/v2);
%         r_hat = x_post;
%         tau2 = v_ext;


        %  method 2 (31)
        
        W_hat = v2*A' / ( v2*(A*A')+sigma2*eye(M) );
        % W_hat = v2*A' / ( v2*diag(diag((A*A')))+sigma2*eye(M) );
        % W_hat = A'/(A*A');
        % W_hat = A';
        W_t = N/trace(W_hat*A)*W_hat;  % Decorrelated
        % W_t = W_t_P(:,:,p);
        H_post = H_hat + W_t*(Y-A*H_hat);
        % B_t = eye(N)-W_t*A;

        % tau2 = 1/N*trace(B_t*B_t')*v2 + 1/N*trace(W_t*W_t')*sigma2; % (31)
        tau2 = 1/N*real(N-2*trace(W_t*A)+trace((W_t'*W_t)*(A*A')))*v2+1/N*trace(W_t'*W_t)*sigma2;
 
        % tau2 = TrBt(p)*v2 + TrWt(p)*sigma2; % (31)
        S_hat = H_post;


        %%  NLE: div-free estimator

        c_mu = vx2*tau2/(tau2 + vx2)* (mux/vx2 + S_hat/tau2);
        C_var =  vx2*tau2/(vx2+tau2);
        L_cal = log(tau2/(tau2+vx2)) + abs(S_hat).^2/tau2 - abs(S_hat-mux).^2/(tau2+vx2);
        pai = rho./( (1-rho).*exp(-L_cal) + rho);
        % pai = repmat(mean(pai,2),[1,T]);  % common sparsity

        eta_mmse = pai.*c_mu;
        var_mmse = pai*C_var + pai.*(1-pai).*abs(c_mu).^2;
        mmseB = mean(real( var_mmse),1 );
        % mmseB =  mean(lambda.*(abs(c_mu).^2+C_var) - abs(eta_mmse).^2); %
        C_t_star = tau2./(tau2-mmseB);
        eta_star = C_t_star.* (eta_mmse - mmseB.*S_hat/tau2);

        H_hat = eta_star;


%         v2  = 1/ (1/mean(mmseB)-1/tau2);  % (39b) % not optimal, also work
%         v2 = max(v2,1e-20);

        v2 = ( norm(Y - A*H_hat,'fro')^2-T*M*sigma2 )/(T*trace(A'*A)); % method-2 (30)
        v2 = max(v2,1e-20);


        H_hat = (1-damp)*H_hat + damp* H_hat_pre;  % add damping
        v2 = (1-damp)*v2 + damp*v2_pre;

        X_hat(:,:,p) = H_hat;
        V2(p) = v2;

        X_hat_pre(:,:,p) = H_hat;    % record last
        V2_pre(p) = v2;


        X_hat_post(:,:,p) = eta_mmse;

        %% EM update
        mux =sum(pai.*c_mu,'all')./sum(pai,'all');  % prior mean
        vx2 =  sum( pai.*(abs(mux-c_mu).^2 + C_var),'all')./sum(pai,'all');  % prior variance
        rho = mean(pai,'all');  % spasity ratio
        if T > 1
            rho = repmat(mean(pai,2),[1,T]); % common sparsity, better for mmv
        end

%         z_mean = A*eta_mmse;
%         z_var = abs(A).^2*var_mmse;
%         z_mean_post = (sigma2*z_mean + Y.*z_var)./(sigma2+z_var);
%         z_var_post = sigma2*z_var./(sigma2+z_var);
%         sigma2 = mean(abs(Y-z_mean_post).^2 + z_var_post,'all');  % better noise estimator
%         % sigma2 = mean(mean(abs(Y-z_mean).^2 ./ (1+z_var/sigma2).^2 + sigma2*z_var./(sigma2+z_var)));  %
%         % sigma2 = mean(mean(abs(Y-z_mean).^2 + z_var)); not right

        Sigma2(p) = sigma2 ;
        Mux(p) = mux ;
        Vx2(p) = vx2 ;
        Rho(:,:,p) = rho ;
    end

    %% record result and EM utilization
    Sigma2 = repmat(mean(Sigma2),[P,1]);
    Mux = repmat(mean(Mux),[P,1]);
    Vx2 = repmat(mean(Vx2),[P,1]);
    Rho = repmat(mean(Rho,3),[1,1,P]);
    
    lambda = mean(mean(Rho,3),2);  % belief indicator
    X_est_mu(:,:,:,it) = X_hat_post; % posterior estimation
           
    MSE_simulation(it) = MSE_simulation(it) + norm(X_hat_post(:) - X(:),2)^2/norm(X(:),2)^2;
    NMSE_iter = norm(X_hat_post(:)-X_hat_pre(:),2)^2/norm(X_hat_pre(:),2)^2;
    if NMSE_iter < 1e-5
        N_iteration = it;
        break;
    else
        N_iteration = niter;
    end
end


end