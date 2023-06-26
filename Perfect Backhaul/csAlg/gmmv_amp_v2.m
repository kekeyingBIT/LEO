function [X_final, lambda_final] = gmmv_amp_v2(Y, Phi, damp, niter, VarN, nns_sel,H_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMMV-AMP algorithm for GMMV CS problem (estimating 3D matrix)
% where incremental EM algorithm is used to learn unknown hyper-parameters.
% An extended version of AMP-NNSPL algorithm in the paper:
% X. Meng et al, "Approximate message passing with nearest neighbor sparsity pattern learning".

% Inputs:
%   Y: received signal
%   Phi: measurement matrix
%   niter: number of AMP iterations
%   tol: termination threshold
%   nns_sel: select nearest neighbor set for NNSPL 0: strctured sparsity   1: clustered sparsity

% Outputs:
%   Xhat: the estimated matrix
%   lambda: belief indicators
%   iter: number of iteration

% Written by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: v2-2019.12.03
% 超参数计算按矢量计算。当目标信号均值方差不同时具有更好的鲁棒性；当均值方差相同时对性能影响待验证
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_CS_model

%% Initializations
[M,Q,P] = size(Y);
[~,N,~] = size(Phi);

% sparsity ratio
alpha = M/N;
normal_cdf = @(x) 1/2.*(1+erf(x/sqrt(2)));
normal_pdf = @(x) 1/sqrt(2*pi).*exp(-x.^2/2);
alpha_grid = linspace(0,10,1024);
rho_SE = (1 - 2/alpha.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)))...
         ./ (1 + alpha_grid.^2 - 2.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)));
lambda = alpha*max(rho_SE).*ones(N,Q,P); 

% noise variance
% snr0 = 100; 
% nvar = sum(abs(Y).^2, 1)/(snr0+1)/M;
% nvar = nvar(ones(1,M),:,:);
nvar = VarN*ones(M,Q,P);

% mean and variance of target signal
xmean = zeros(N,Q,P);
Phi_fro = sum(sum(abs(Phi).^2, 1), 2);
Phi_fro = Phi_fro(:, ones(1,Q), :);
xvar = (sum(abs(Y).^2, 1) - M.*nvar(1,:,:)) ./ Phi_fro ./ (alpha*max(rho_SE));
xvar = xvar(ones(1,N),:,:);

% other parameters initialization
Xhat = xmean; 
v = xvar;
V = ones(M,Q,P);
Z = Y;
% nearest neighboor sparsity pattern learning
index = ones(N,Q,P);
index_left = [zeros(N,1,P), index(:,1:Q-1,:)];
index_right = [index(:,2:Q,:), zeros(N,1,P)];
index_ahead = cat(3, zeros(N,Q,1), index(:,:,1:P-1));
index_latter = cat(3, index(:,:,2:P), zeros(N,Q,1));
index_3D = index_left + index_right + index_ahead + index_latter; % index_up + index_down + 
clear index
clear index_left
clear index_right
clear index_ahead
clear index_latter
% allocate memory
D = zeros(N,Q,P);
C = zeros(N,Q,P);
L_cal = zeros(N,Q,P);
pai = zeros(N,Q,P);
A = zeros(N,Q,P);
B = zeros(N,Q,P);
Yhat = zeros(M,Q,P);

NMSE = zeros(1,niter);
NMSE_re = zeros(1,niter);
Xhat_pre = zeros(N,Q,P);


tol = 1e-5;
tol_re = 1e-4;
%% AMP iteration
for iter = 1:niter
    %Xhat_pre = Xhat;
    V_pre = V;
    for p = 1:P
        % factor node update 
        V(:,:,p) = damp.*V_pre(:,:,p) + (1-damp).*abs(Phi(:,:,p)).^2*v(:,:,p);
        Z(:,:,p) = damp.*Z(:,:,p) + (1-damp).*...
                   (Phi(:,:,p)*Xhat(:,:,p) - (Y(:,:,p)-Z(:,:,p))./(nvar(:,:,p)+V_pre(:,:,p)).*V(:,:,p));
                                          
        % variable node update 
        D(:,:,p) = 1 ./ ((abs(Phi(:,:,p)).^2).'*(1./(nvar(:,:,p)+V(:,:,p))));
        C(:,:,p) = Xhat(:,:,p) + D(:,:,p).*(Phi(:,:,p)'*((Y(:,:,p)-Z(:,:,p))./(nvar(:,:,p)+V(:,:,p))));
        
        % compute posterior mean and variance
        L_cal(:,:,p) = (1/2).*(log(D(:,:,p)./(D(:,:,p)+xvar(:,:,p))) + abs(C(:,:,p)).^2./D(:,:,p) - ...
                                   abs(C(:,:,p)-xmean(:,:,p)).^2./(D(:,:,p)+xvar(:,:,p))); 
        pai(:,:,p) = lambda(:,:,p) ./ (lambda(:,:,p)+(1-lambda(:,:,p)).*exp(-L_cal(:,:,p))); 
        A(:,:,p) = (xvar(:,:,p).*C(:,:,p)+xmean(:,:,p).*D(:,:,p)) ./ (D(:,:,p)+xvar(:,:,p));
        B(:,:,p) = (xvar(:,:,p).*D(:,:,p)) ./ (xvar(:,:,p)+D(:,:,p));
        Xhat(:,:,p) = pai(:,:,p).*A(:,:,p);
        v(:,:,p) = pai(:,:,p).*(abs(A(:,:,p)).^2+B(:,:,p)) - abs(Xhat(:,:,p)).^2;
        
        % reconstruct received signal
        Yhat(:,:,p) = Phi(:,:,p)*Xhat(:,:,p);
    end
    
    % hyper-parameter update based on EM learning
    xmean = sum(pai.*A, 1)./sum(pai, 1);
    xmean = xmean(ones(1,N),:,:);
    xvar = sum(pai.*(abs(xmean-A).^2+B), 1)./sum(pai, 1); 
    xvar = xvar(ones(1,N),:,:);
%     nvar = sum(abs(Y-Z).^2./abs(1+V./nvar).^2 + V./(1+V./nvar))/M;
%     nvar = nvar(ones(1,M),:,:);
    
    % refine the update rule for sparsity ratio
    switch nns_sel
        case 0 % structured sparsity(common support)
            pai_tmp = sum(sum(pai,3),2)./Q./P;
            lambda = pai_tmp(:, ones(1,Q), ones(1,P));
        case 1  % clustered sparsity
            pai_left = [zeros(N,1,P), pai(:,1:Q-1,:)];
            pai_right = [pai(:,2:Q,:), zeros(N,1,P)];
            pai_ahead = cat(3, zeros(N,Q,1), pai(:,:,1:P-1));
            pai_latter = cat(3, pai(:,:,2:P), zeros(N,Q,1));
            lambda = (pai_left + pai_right + pai_ahead + pai_latter) ./ index_3D;

    end
    
    % Xhat_store(:,:,:,iter) = Xhat;
    % stopping critrria
    NMSE_iter =  norm(Xhat(:)-Xhat_pre(:))^2 / norm(Xhat_pre(:))^2;
    NMSE_re_iter = norm(Y(:)-Yhat(:))^2 / norm(Y(:))^2;
    NMSE(iter) = (norm(Xhat(:)-H_true(:),2)/norm(H_true(:),2))^2;
    NMSE_re(iter) = (norm(Yhat(:)-Y(:),2) / norm(Y(:),2))^2;
    if NMSE_iter < tol || NMSE_re_iter < tol_re
        X_final =  Xhat;
        lambda_final = lambda;
        break;
    else
        X_final =  Xhat;
        lambda_final = lambda;
    end
    Xhat_pre = Xhat;

         
end
% figure
% plot(10*log10(NMSE),'-.');hold on
% [val,idx] = min(10*log10(NMSE(1:iter)));
% X_final =  Xhat_store(:,:,:,idx);
end