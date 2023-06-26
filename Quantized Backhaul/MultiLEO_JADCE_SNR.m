function [H_hat_final, est_idx_final] = MultiLEO_JADCE_SNR(ce_choice,ob_nIBI,Phi,Nvar,Agc, ...
    H_true,Act_true,ch,leo,user,ob,snr,ch_para)
% perform multi-satellite and multi-block JADCE

% system parameter
N_u = user.N_u;
L_max = ch.L_max;
N_r = prod(leo.N_r);
N_p = leo.N_p;
N_block = size(ob_nIBI,3);

% % algorithm hyper-parameter 
% th_cg = 0.5;  % set high to  avoid false alarm
% epsilon = 0.02; 
% 
% iterT = 40;
% damp = 0.3;

% algorithm hyper-parameter 
th_cg = 0.5;  % set high to  avoid false alarm, set low to to avoid miss detection
epsilon = 0.02; 
if ob == 1
    iterT = 50;
    damp = 0.5;
elseif ob == 2 
    iterT = 40;
    damp = 0.4;
else
    iterT = 40;
    damp = 0.2;
end

struct = 0;
% JADCE  main body

switch ce_choice  % select signal processing method
    case 'OAMP-0'
        H_hat_all = zeros(N_u*L_max,N_r,N_p);
        est_idx_all = zeros(N_u,N_p);
        for np = 1:N_p
            agc = Agc(np);
            % agc =  1;
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            H_mtx = agc* repmat(H_true(:,:,np),[1,1,N_block]); % assume each block experience same channel->also can be different
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            % [H_est, lambda] = gmmv_oamp(ob_nIBI_np,Phi,damp, iterT,nvar,H_mtx);
            [H_est, lambda,N_iteration] = gmmv_oamp_app(ob_nIBI_np,Phi,damp, iterT,nvar,H_mtx);
           
            H_hat =  mean(H_est(:,:,:,N_iteration)/agc,3); % if same channel is used for each block

            epsilon_cg = epsilon*max(abs(H_hat(:)));

            supp = ( sum(abs(H_hat)>epsilon_cg,2)/N_r ) > th_cg;
            est_idx = double( sum( reshape(supp,[L_max,N_u]), 1) > 0 ).';
            H_hat_all(:,:,np) = H_hat;
            est_idx_all(:,np) = est_idx;
        end
        H_hat_final = H_hat_all;
        est_idx_final = double(mean(est_idx_all,2) > 0); % if one LEO detect active, it is active, (easy to result high false alarm)
        
    case 'OAMP-1' % simplified complexity
        H_hat_all = zeros(N_u*L_max,N_r,N_p);
        est_idx_all = zeros(N_u,N_p);
        [M,N,P] = size(Phi);
        U_all = zeros(M,M,P);
        S_all = zeros(M,N,P);
        V_all = zeros(N,N,P);
        for p = 1:P
            A = Phi(:,:,p);
            [U0,S0,V0] = svd(A); % A =  U*\Sigma*U'
            U_all(:,:,p) = U0;
            S_all(:,:,p) = S0;
            V_all(:,:,p) = V0;
        end
        
        for np = 1:N_p
            agc = Agc(np);
            % agc =  1;
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            H_mtx = agc* repmat(H_true(:,:,np),[1,1,N_block]); % assume each block experience same channel->also can be different
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            % [H_est, lambda] = gmmv_oamp(ob_nIBI_np,Phi,damp, iterT,nvar,H_mtx);
            [H_est, lambda,N_iteration] = gmmv_oamp_app2(ob_nIBI_np,Phi,damp, iterT,nvar,H_mtx,...
                U_all,S_all,V_all);
           
            H_hat =  mean(H_est(:,:,:,N_iteration)/agc,3); % if same channel is used for each block

            epsilon_cg = epsilon*max(abs(H_hat(:)));

            supp = ( sum(abs(H_hat)>epsilon_cg,2)/N_r ) > th_cg;
            est_idx = double( sum( reshape(supp,[L_max,N_u]), 1) > 0 ).';
            H_hat_all(:,:,np) = H_hat;
            est_idx_all(:,np) = est_idx;
        end
        H_hat_final = H_hat_all;
        est_idx_final = double(mean(est_idx_all,2) > 0); % if one LEO detect active, it is active, (easy to result high false alarm)
    case 'AMP'
        H_hat_all = zeros(N_u*L_max,N_r,N_p);
        est_idx_all = zeros(N_u,N_p);
        for np = 1:N_p
            agc = Agc(np);
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            H_mtx = agc* repmat(H_true(:,:,np),[1,1,N_block]); % assume each block experience same channel->also can be different
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            [H_est,lambda]=gmmv_amp_v2(ob_nIBI_np, Phi, damp, iterT,nvar,struct,H_mtx);  % assume same pilot is adopt for each stream

            H_hat =  mean(H_est/agc,3); % if same channel is used for each block

            epsilon_cg = epsilon*max(abs(H_hat(:)));

            supp = ( sum(abs(H_hat)>epsilon_cg,2)/N_r ) > th_cg;
            est_idx = double( sum( reshape(supp,[L_max,N_u]), 1) > 0 ).';
            H_hat_all(:,:,np) = H_hat;
            est_idx_all(:,np) = est_idx;
        end
        H_hat_final = H_hat_all;
        est_idx_final = double(mean(est_idx_all,2) > 0); % if one LEO detect active, it is active, (easy to result high false alarm)

%         % debug verify
%         nmse = ( norm(H_hat_final(:)-H_true(:),2) / norm(H_true(:),2) )^2;
%         acc = ( sum( xor( est_idx_final,Act_true ) ) == 0 );

    case 'OAMP-2'  % each leo perform JADCE independently
         H_hat_all = zeros(N_u*L_max,N_r,N_p);
         H_hat_all_refined = zeros(N_u*L_max,N_r,N_p);
         est_idx_all = zeros(N_u,N_p);
         [M,N,P] = size(Phi);
         U_all = zeros(M,M,P);
         S_all = zeros(M,N,P);
         V_all = zeros(N,N,P);
         for p = 1:P
             A = Phi(:,:,p);
             [U0,S0,V0] = svd(A); % A =  U*\Sigma*U'
             U_all(:,:,p) = U0;
             S_all(:,:,p) = S0;
             V_all(:,:,p) = V0;
         end
         
         for np = 1:N_p
            agc = Agc(np);
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            H_mtx = agc* repmat(H_true(:,:,np),[1,1,N_block]); % assume each block experience same channel->also can be different
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            % [H_est, lambda] = gmmv_oamp(ob_nIBI_np,Phi,damp, iterT,nvar,H_mtx);
            [H_est, lambda,N_iteration] = gmmv_oamp_app2(ob_nIBI_np,Phi,damp, iterT,nvar,H_mtx,...
                U_all,S_all,V_all);
            
            H_hat =  mean(H_est(:,:,:,N_iteration),3); % if same channel is used for each block
            
            epsilon_cg = epsilon*max(abs(H_hat(:)));
            supp = ( sum(abs(H_hat)>epsilon_cg,2)/N_r ) > th_cg; 
            est_idx = double( sum( reshape(supp,[L_max,N_u]), 1) > 0 ).';

            % DoA = EstDoA(H_hat,supp,est_idx,leo.N_r,ch.L_max);
            DoA = EstDoA2(H_hat,supp,est_idx,leo.N_r,ch.L_max);
            [H_hat_refined,es] = chRebuild(H_hat,DoA,est_idx,supp,ch,leo,user,ch_para, H_true);
            % [H_hat_refined,es] = chRebuild(agc*H_true(:,:,np),DoA,est_idx,supp,ch,leo,user,ch_para, H_true);
            H_hat_refined = H_hat_refined/agc;

            H_hat_all(:,:,np) = H_hat;
            est_idx_all(:,np) = est_idx;
            H_hat_all_refined(:,:,np) = H_hat_refined;

         end
        
         

         H_hat_final = H_hat_all_refined;
         est_idx_final = double(mean(est_idx_all,2) > 0); % if one LEO detect active, it is active, (easy to result high false alarm)
        
%          % debug verify
%          nmse1 = 10*log10( ( norm(H_hat_all(:)-H_true(:),2) / norm(H_true(:),2) )^2)
%          nmse2 = 10*log10(( norm(H_hat_all_refined(:)-H_true(:),2) / norm(H_true(:),2) )^2)
%          acc = ( sum( xor( est_idx_final,Act_true ) ) == 0 )

    case 'Proposed-2' % utilize neighboring infomation -> no gain
         H_hat_all = zeros(N_u*L_max,N_r,N_p);
         H_hat_all_refined = zeros(N_u*L_max,N_r,N_p);
         est_idx_all = zeros(N_u,N_p);
         supp_all =  zeros(N_u*L_max,N_p);

         H_hat_final = zeros(N_u*L_max,N_r,N_p);
         % independent channel estimation at each LEO
        for np = 1:N_p
            agc = Agc(np);
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            H_mtx = agc* repmat(H_true(:,:,np),[1,1,N_block]); % assume each block experience same channel->also can be different
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            [H_est,lambda]=gmmv_amp_v2(ob_nIBI_np, Phi, damp, iterT,struct,H_mtx);  % assume same pilot is adopt for each stream
            
            H_hat =  mean(H_est/agc,3); % if same channel is used for each block
            
            epsilon_cg = epsilon*max(abs(H_hat(:)));
            supp = ( sum(abs(H_hat)>epsilon_cg,2)/N_r ) > th_cg; 
            est_idx = double( sum( reshape(supp,[L_max,N_u]), 1) > 0 ).';

            DoA = EstDoA(H_hat,supp,est_idx,leo.N_r,ch.L_max);
            theta_true = ch_para.theta(np,:).';
            phi_true = ch_para.phi(np,:).';
            est_set = find(est_idx>0);
            act_set = find(Act_true>0);
            [H_hat_refined,es] = chRebuild(H_hat,DoA,est_idx,supp,ch,leo,user,ch_para, H_true);
            H_hat_all(:,:,np) = H_hat;
            est_idx_all(:,np) = est_idx;
            H_hat_all_refined(:,:,np) = H_hat_refined;
%             supp_all(:,np) = supp;
% 
%             nmse1 = ( norm(H_hat-H_true(:,:,np),'fro') / norm(H_true(:,:,np),'fro') )^2
%             nmse2 = ( norm(H_hat_refined-H_true(:,:,np),'fro') / norm(H_true(:,:,np),'fro') )^2
         end
         
         nmse1 = 10*log10( ( norm(H_hat_all(:)-H_true(:),2) / norm(H_true(:),2) )^2);
         nmse2 = 10*log10(( norm(H_hat_all_refined(:)-H_true(:),2) / norm(H_true(:),2) )^2);

         % merged side-aud information from  neighboring satellite
         for np = 1:N_p
             est_idx_tp = est_idx_all(:,np); % each AP's own detected information
             est_set_tp = find(est_idx_tp>0);
             add_set_tp = find( sum(est_idx_all,2) > 0 ) ; % 
             miss_idx_tp = setdiff(add_set_tp, est_set_tp);    % miss alarm set
             H_tp = H_hat_all_refined(:,:,np);

            
             % SIC for enhanced estimation for miss-detected user
             
             supp_res_mat =( (miss_idx_tp-1)*ch.L_max + repmat([1:ch.L_max],[length(miss_idx_tp),1]) ).';
             supp_res = supp_res_mat(:);

             H_est_ave = H_tp(supp_res,:);
             H17 = H_true(supp_res,:,np);
            
             % searching ESPRIT for miss detected user -> inaccurate method
             % due to too few snapshot
             R_res = diag(H_est_ave*H_est_ave');
             epsilon_res = epsilon*max(abs(R_res));
             supp_res_act = R_res>epsilon_res;
             est_res_act = double( sum( reshape(supp_res_act,[L_max,length(miss_idx_tp)]), 1) > 0 ).';
             DoA_res = EstDoA(H_est_ave, supp_res_act,est_res_act,leo.N_r,ch.L_max);
             [H_res_refined,es] = chRebuild(H_est_ave,DoA_res,est_res_act,supp_res_act,ch,leo,user,ch_para, H_true);

             H_tp(supp_res,:) = H_res_refined;
             H_tp_true = H_true(supp_res,:,np);
             n1 = ( norm(H_est_ave(:)-H_tp_true(:),2) / norm(H_tp_true(:),2) )^2;
             n2 = ( norm(H_res_refined(:)-H_tp_true(:),2) / norm(H_tp_true(:),2) )^2;

             H_hat_final(:,:,np) = H_tp;
             % est_idx_final

             Pe(np) = sum( xor( est_idx_tp,Act_true ) )/length(Act_true);
             
         end    
                  
         nmse1 = 10*log10( ( norm(H_hat_all(:)-H_true(:),2) / norm(H_true(:),2) )^2);
         nmse2 = 10*log10(( norm(H_hat_all_refined(:)-H_true(:),2) / norm(H_true(:),2) )^2);
         nmse3 = 10*log10(( norm(H_hat_final(:)-H_true(:),2) / norm(H_true(:),2) )^2);

         % output result
         H_hat_final = H_hat_all;
         est_idx_final = double(mean(est_idx_all,2) > 0); % if one LEO detect active, it is active, (easy to result high false alarm)        
         Pe_refined(np) = sum( xor( est_idx_final,Act_true ) )/length(Act_true);
         acc = ( sum( xor( est_idx_final,Act_true ) ) == 0 );
         
    case 'Perfect-CIR'
        H_hat_final = H_true;
        est_idx_final = Act_true;


    case 'SOMP'    % need to be correct
        H_hat_all = zeros(N_u*L_max,N_r,N_p);
        est_idx_all = zeros(N_u,N_p);
        for np = 1:N_p
            agc = Agc(np);
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            niterN = length(find(Act_true>0))*ch.N_c;
            [H_est, supp_set, ~] = gmmv_omp(ob_nIBI_np, Phi, nvar,niterN);
            H_hat =  mean(H_est/agc,3); % if same channel is used for each block
            
            supp = zeros(L_max*N_u,1);
            supp(supp_set,:) = ones(length(supp_set),1);
            est_idx = double( sum( reshape(supp,[L_max,N_u]), 1) > 0 ).';
            H_hat_all(:,:,np) = H_hat;
            est_idx_all(:,np) = est_idx;
        end
        H_hat_final = H_hat_all;
        est_idx_final = double(mean(est_idx_all,2) > 0); % if one LEO detect active, it is active, (easy to result high false alarm)

        nmse = ( norm(H_hat_final(:)-H_true(:),2) / norm(H_true(:),2) )^2;
        acc = ( sum( xor( est_idx_final,Act_true ) ) == 0 );

    case {'Oracle-LS','Oracle-LMMSE'}  % need to be correct

        H_hat_all = zeros(N_u*L_max,N_r,N_p);
        for np = 1:N_p
            agc = Agc(np);
            ob_nIBI_np = agc* ob_nIBI(:,:,:,np);
            nvar = agc^2 * Nvar(np); % equivalent noise variance after agc
            supp =  find(sum(abs(H_true(:,:,np)),2)>0);
            H_supp = H_true(supp,:,np);
            H_est = Oracle_LS_3D(ob_nIBI_np, Phi,supp,nvar,ce_choice,H_supp);
            H_hat =  mean(H_est/agc,3); % if same channel is used for each block
            H_hat_all(:,:,np) = H_hat;

        end
        H_hat_final = H_hat_all;
        est_idx_final = Act_true; % if one LEO detect active, it is active, (easy to result high false alarm)
        
        nmse = ( norm(H_hat_final(:)-H_true(:),2) / norm(H_true(:),2) )^2;
        acc = ( sum( xor( est_idx_final,Act_true ) ) == 0 );


    otherwise
        error('Undefined JADCE algorithm !');


end



end

