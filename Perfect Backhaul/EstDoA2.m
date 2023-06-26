function  DoA = EstDoA2(H_est,supp,est_idx,Nr,Lmax)
% An implementation of "closed-form 2-D Angle Estimation with Rectangular in 
% Element Space or Beamspace via Unitary ESPRIT"，1996

% v2, add spatial smooth
%% paramaters setting

act_set = find(est_idx==1);
DoA = zeros(length(act_set),2);
N = Nr(1);
M = Nr(2);

remain_x = 3; remain_y = 3;
Mx_sub = N-remain_x; Gx = N-Mx_sub+1;
My_sub = M-remain_y; Gy = M-My_sub+1;
%% special matrix define

% define PI matrix
   % PI = @(K) flip(eye(K));
% define Q matrix
    Q_N  = ConQ(Mx_sub); Q_M  =ConQ(My_sub);
% define J matrix
%     J1  = [eye(N-1),zeros(N-1,1)];
    J2_N =  [zeros(Mx_sub-1,1),eye(Mx_sub-1)];
    J2_M =  [zeros(My_sub-1,1),eye(My_sub-1)];
% define K matrix
    tmp = (ConQ(Mx_sub-1))'*J2_N*ConQ(Mx_sub);
    K1 = real(tmp);
    K2 = imag(tmp);
    tmp2 = (ConQ(My_sub-1))'*J2_M*ConQ(My_sub);
    K3 = real(tmp2);
    K4 = imag(tmp2);
% 2D - ESPRIT 

    for k = 1:length(act_set)
        % Signal transmission
        u = act_set(k);
        seq_idx = (u-1)*Lmax+1:u*Lmax;
        % supp_u  = find ( supp(seq_idx) > 0)   ;
        Hk = H_est(seq_idx,: );

        % X_n = (Hk(supp_u,:)).';
        
        X_n = (Hk(supp(seq_idx) > 0,:)).';
        snp = size(X_n,2);
        
        % spatial smoothing
        X_n_tilde = zeros(Mx_sub*My_sub,Gx*Gy*snp );
        
        for gy = 1:Gy
            for gx = 1:Gx
                % J_xy = kron(SelMat(gx,Gx,Mx_sub), SelMat(gy,Gy,My_sub));
                J_xy = kron(SelMat(gy,Gy,My_sub),SelMat(gx,Gx,Mx_sub));
                col_idx = (gy-1)*Gx*snp+ [(gx-1)*snp+1: gx*snp];
                X_n_tilde(:,col_idx) = J_xy*X_n;
            end
        end
        % 2D u-esprit
        % step 1
        Y  = kron(Q_M',Q_N')*X_n_tilde;
        % Y  = kron(Q_M',Q_N')*X_n;
        Z = [real(Y),imag(Y)];
        [U,~,~] = svd(Z);
        Es = U(:,1:1);
        % step 2
        K_mu1 = kron(eye(My_sub),K1);
        K_mu2 = kron(eye(My_sub),K2);
        K_v1 = kron(K3,eye(Mx_sub));
        K_v2 = kron(K4,eye(Mx_sub));

        Psi_mu =(K_mu1*Es)\K_mu2*Es;
        Psi_v =(K_v1*Es)\K_v2*Es;
        Psi = Psi_mu+1i*Psi_v;
        Omi = eig(Psi);
        mux_hat  = 2*atan(real(Omi));
        muy_hat = 2*atan(imag(Omi));
        theta_hat = atan2(muy_hat,mux_hat) -pi;
        phi_hat = asin(1/pi*sqrt((mux_hat.^2+muy_hat.^2)));
        DoA(k,1) = phi_hat .* 180./pi;
        % DoA(k,2) = theta_hat .* 180./pi ;
        DoA(k,2) = theta_hat .* 180./pi;
    end

end  % function end


function Q_N = ConQ(N)
    PI = @(K) flip(eye(K));
    if mod(N,2) == 0
    K  = N/2;
    Q_N = 1/sqrt(2)*[eye(K), 1i*eye(K);PI(K),-1i*PI(K)];
    else
    K = (N-1)/2;
    Q_N = 1/sqrt(2)*...
        [eye(K), zeros(K,1), 1i*eye(K);...
        zeros(1,K),sqrt(2),zeros(1,K);...
        PI(K),zeros(K,1),-1i*PI(K)];
    end
end

function J_s = SelMat(gs,Gs,Ms_sub)
% selection matrix for spatial smoothing
J_s = [zeros(Ms_sub,gs-1), eye(Ms_sub),zeros(Ms_sub,Gs-gs)];

end