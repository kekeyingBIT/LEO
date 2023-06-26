function [rec_pilot] = calPilotIBI(TS_symbol,Hd_alg_est,act_alg_est,N_block,fr)

[Ndim,N_r,N_p] = size(Hd_alg_est);
N_u = size(act_alg_est,1);
L_max = Ndim/N_u;

n_symbol_pblock = fr.n_symbol_pblock;
n_ts = fr.n_ts_pblock;
n_all_pblock = fr.n_all_pblock;

est_set = find(act_alg_est>0);
N_a = length(est_set);
% N_u = size(TS_symbol,2);


rec_pilot = zeros(n_all_pblock*N_block+L_max-1,N_r,N_p);

%% TS seq generated at receier locally
for u = 1:N_a
    user_u = est_set(u);
    Data_tmp_sym_t = zeros(n_symbol_pblock,1);
    
    PilotFrame = zeros(n_all_pblock*N_block,1);
    for blk = 1:N_block
        TS_tmp_sym_t = TS_symbol((blk-1)*n_ts+1:blk*n_ts,user_u);
        sig_u_blk = [TS_tmp_sym_t; Data_tmp_sym_t];
        PilotFrame((blk-1)*n_all_pblock+1:blk*n_all_pblock,:) = sig_u_blk;
    end  % block end

    for np = 1:N_p
        % generate rec TS
        h_u = Hd_alg_est((user_u-1)*L_max+1:user_u*L_max,:,np);

        for nr = 1:N_r
            rec_pilot(:,nr,np) = rec_pilot(:,nr,np) + conv(PilotFrame,h_u(:,nr));
        end  % N_r end

    end % N_p end

end  % N_a end

end