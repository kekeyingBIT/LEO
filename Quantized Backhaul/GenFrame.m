function [sframe_all,sframe_symbol,mask_pos_all,fr] = GenFrame(act_set,TS_symbol,Data_bits,N_block,M_ord,mask,fq)

mask_cr = mask.cr;  
mask_set = mask.set;
N_bits = log2(M_ord); % bits number per symbol
N_a = length(act_set);
N_k =  fq.N_k;

% data frame without considering mask

% frame level parameter
[ts_len_pframe, N_u] = size(TS_symbol);
symbol_len_pframe = size(Data_bits,1)/N_bits;
N_s = size(Data_bits,2);   % data stream number for each tx

% block level parameter after considering spread 
n_ts_pblock = floor(ts_len_pframe/N_block);
n_symbol_subpblock =floor(symbol_len_pframe/N_block); 
n_symbol_pblock = n_symbol_subpblock/mask_cr;  % spread info length in each block
n_all_pblock =  n_ts_pblock+n_symbol_pblock;  

% allocation memory
sframe_all = zeros(ts_len_pframe + symbol_len_pframe/mask_cr,N_s,N_u);
sframe_symbol = zeros(symbol_len_pframe/mask_cr,N_s,N_u);

% introduce masking to the input symbol
[mask_cw_len,mask_cb_size] = size(mask_set);
mask_nz_num = length(mask_cw_len*mask_cr);
mask_blk_num = n_symbol_subpblock/mask_cw_len;
mask_pos_all = zeros(mask_blk_num*mask_nz_num,N_u);
% mask  pattern generation
for uu = 1:N_u
    cw_sel = randi(mask_cb_size,1,mask_blk_num); % different mask block adopt different mask
    mask_sc_pos_tmp = mask_set(:,cw_sel);
    mask_sc_pos_idx = mask_sc_pos_tmp(:);
    mask_sc_pos = find(mask_sc_pos_idx == 1);
    mask_pos_all(:,uu) = mask_sc_pos;
end




% DFT-s-OFDM ,mapping
for u = 1:N_a
    act_idx = act_set(u);
    frame_all_u= zeros(N_block*n_all_pblock,N_s);
    frame_symbol_u = zeros(N_block*n_symbol_pblock,N_s);
    TS_tmp_u = TS_symbol(:,act_idx);
    dataEnc_u = Data_bits(:,:,act_idx);    
    mask_sc_pos = mask_pos_all(:,act_idx);
    
    for s = 1:N_s  % different data stream
        dataEnc_u_s = dataEnc_u(:,s);
        for i = 1:N_block
            dataEnc_u_s_info = dataEnc_u_s((i-1)*n_symbol_subpblock*N_bits+1:i*n_symbol_subpblock*N_bits);
            TS_tmp_sym_t = TS_tmp_u((i-1)*n_ts_pblock+1:i*n_ts_pblock);
            n_data_psym_sc = n_symbol_subpblock*mask_cr;  % data len after spead coding /masking
            Fo_DFT = sqrt(1/n_data_psym_sc)*dftmtx(n_data_psym_sc);
            
            Data_tmp_sym_t = zeros(N_k/mask_cr,1);
            for sblk = 1:double(1/mask_cr)
                Data_sym_f_info = dataEnc_u_s_info((sblk-1)*n_data_psym_sc*N_bits+1:sblk*n_data_psym_sc*N_bits,:);

                Data_symbol_mod_info = qammod(Data_sym_f_info,M_ord,'gray','InputType', 'bit','UnitAveragePower', true);

                % DFT-s-OFDM predure
                % DFT-operation (K-point DFT)
                Data_symbol_mod_info = Fo_DFT'*Data_symbol_mod_info;
                % spread-operation (Subcarrier allocation)
                Data_symbol_mod = zeros(N_k,1);
                Data_symbol_mod(mask_sc_pos,:) = Data_symbol_mod_info;
                % OFDM operation (M-point IFFT)
                Data_tmp_sym_t_sblk = 1/sqrt(N_k)*dftmtx(N_k)*Data_symbol_mod;
                Data_tmp_sym_t ((sblk-1)*N_k+1:sblk*N_k,:) = Data_tmp_sym_t_sblk;
                
            end
            
            Block_symbol = [TS_tmp_sym_t; Data_tmp_sym_t];  % assembled symbol in one block
            frame_all_u((i-1)*n_all_pblock+1:i*n_all_pblock,s) =  Block_symbol;

            % In TD/FD, this variable stands for different meaning
            frame_symbol_u((i-1)*n_symbol_pblock+1:i*n_symbol_pblock,s) = Data_tmp_sym_t;

        end  % block end
    end % Ns end

    sframe_all(:,:,act_idx) = frame_all_u;
    sframe_symbol(:,:,act_idx) = frame_symbol_u;
end % N_a end

% record block parameter for use
fr.n_ts_pblock = n_ts_pblock;
fr.n_symbol_subpblock = n_symbol_subpblock;
fr.num_subblock = 1/mask_cr;
fr.n_symbol_pblock = n_symbol_pblock;
fr.n_all_pblock = n_all_pblock;

end % function end

