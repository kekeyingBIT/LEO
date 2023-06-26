function [Ht_all_fre] = GenFreqCh(Hd_alg_est,act_alg_est,mask_pos_all,fr)
% generate frequency domain channel



%% parameters extraction
[Ndim,N_r,N_p] = size(Hd_alg_est);
N_u = size(act_alg_est,1);
L_max = Ndim/N_u;

n_symbol_pblock = fr.n_symbol_pblock;
n_ts = fr.n_ts_pblock;
n_all_pblock = fr.n_all_pblock;
n_symbol_subpblock = fr.n_symbol_subpblock;
num_subblock = fr.num_subblock;

est_set = find(act_alg_est>0);
N_a = length(est_set);

N_s= n_symbol_subpblock;

Fo = 1/sqrt(N_s)*dftmtx(N_s);

%% build the detection channel matrix

% frequency domain equalization by DFT

Ht_all_fre = zeros(N_r,N_a,N_s,N_p);
% Ht_all_fre_mat = zeros(N_r*N_s,N_act*N_s);
% tic
for np = 1:N_p
    for u = 1:N_a
        user_u = est_set(u);
        user_seq = (user_u-1)*L_max+1:user_u*L_max;
        hu = Hd_alg_est(user_seq,:,np);
        for nr = 1:N_r
            hu_nr = (hu(:,nr)).';
            % hrow = [hu_nr(1), zeros(1,N_s-L_max), flip(hu_nr(2:end))];
            % ht_sc_u = Fo'*toeplitz([hrow(1) fliplr(hrow(2:end))],hrow)*Fo;
            % fre_sc_u   = diag(ht_sc_u);
            
            % delete the blank channel
            mask_pos_u = mask_pos_all(:,user_u);
            mask_zero_pos_u = setdiff((1:N_s),mask_pos_u);
            
            fre_sc_u = conj(fft(hu_nr',N_s));
            fre_sc_u(mask_zero_pos_u) = zeros(length(mask_zero_pos_u),1);
            Ht_all_fre(nr,u,:,np) = fre_sc_u;

        end
    end
end
% toc


end

