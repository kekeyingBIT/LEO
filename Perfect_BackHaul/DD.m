function [data_recv_info_bit,data_recv_info_bit_Uncoded] = DD(rec_sig,PilotIBI,...
    chMtx1,est_idx,mask_pos_all,dd_sel,para,Coded_Sel,...
    Coded_Info_bits,Uncoded_Info_bits,sframe_symbol)
%% parse parameter 
est_set = find(est_idx > 0);
[N_r,N_a,N_k,N_p] = size(chMtx1);  % []
N_u = length(est_idx);

LDPC = para.LDPC;
mask = para.mask;
N_block = para.N_block;
L_max = para.L_max;
M_ord  = para.M_ord;
Nvar = para.Nvar;
N_bits = log2(M_ord); 

fr = para.fr;
n_symbol_pblock = fr.n_symbol_pblock;
n_ts = fr.n_ts_pblock;
n_symbol_subpblock = fr.n_symbol_subpblock;
N_subblock = fr.num_subblock;
n_all_pblock = fr.n_all_pblock;

% DFT matrix predefine
N_k = n_symbol_subpblock;
N_b = n_symbol_pblock; 

N_info = N_k*mask.cr;
Fo_DFT = sqrt(1/N_info)*dftmtx(N_info);  % IDFT length of DFT-s-OFDM
Fo_FFT = 1/sqrt(N_k)*dftmtx(N_k);   % Subcarrier number

% memory allocate  
N_info_pblock_bit = N_subblock*N_info*N_bits;
data_recv_info = zeros(N_info_pblock_bit,N_block,N_a,N_p);
switch Coded_Sel
    case 'LDPC'
        data_recv_info_bit = zeros(LDPC.numInfBits,N_u);
        data_recv_info_bit_Uncoded = zeros(LDPC.numTotBits,N_u);
    case 'No'
        data_recv_info_bit = zeros(LDPC.numTotBits,N_u);
        data_recv_info_bit_Uncoded = zeros(LDPC.numTotBits,N_u);
end
% nzp_ = 

%% data detection
for np = 1:N_p
    for blk =  1:N_block
        %% interference cancellation part
        ts_idx = (blk-1)*N_b + ((blk-1)*n_ts+1:blk*n_ts);
        data_idx = blk*n_ts + ( (blk-1)*N_b +1 :blk*N_b );

        If_bf_idx = ts_idx(end)+1:ts_idx(end)+L_max-1;
        If_af_idx = data_idx(end)+1:data_idx(end) + L_max-1;

        bf_blk = [ PilotIBI(If_bf_idx,:,np); zeros(N_b-(L_max-1),N_r)];
        af_blk = [ PilotIBI(If_af_idx,:,np); zeros(N_b-(L_max-1),N_r)];
        cir_data_blk = [ rec_sig(If_af_idx,:,np); zeros(N_b-(L_max-1),N_r)];

        data_sym = rec_sig(data_idx,:,np);
        cir_data_sym = data_sym  + cir_data_blk- bf_blk- af_blk;

        %  demodulate of the time domain symbol
        x_nc_block = zeros(N_a,N_k, N_subblock);
        for sblk = 1 : N_subblock
            data_subblock = cir_data_sym((sblk-1)*N_k+1:sblk*N_k,:);
            fre_data_sym = Fo_FFT'*data_subblock;  % data in sblk-th subblock

            % prepare the mask matrix and memory
            Ht_all_fre = chMtx1(:,:,:,np);
            x_nc_subblock = zeros(N_a,N_k);
            pos_mat = zeros(N_k,N_a);
            for n_a = 1:N_a
                current_idx = mask_pos_all(:,est_set(n_a));
                pos_mat(current_idx,n_a) = ones(length(current_idx),1); % index matrix regeneration
            end

            for k = 1:N_k
                y_s = fre_data_sym(k,:).';
                channel_s = Ht_all_fre(:,:,k);

                s_u_equ = find(pos_mat(k,:)==1);
                channel_s_equ =  channel_s(:,s_u_equ);

                % LS
                % x_s_u_fre = (channel_s_equ'*channel_s_equ)\channel_s_equ'*y_s;
                % LMMSE
                x_s_u_fre = (channel_s_equ'*channel_s_equ + Nvar(np)*mask.cr*eye(size(channel_s_equ,2)))\channel_s_equ'*y_s;

                x_nc_subblock(s_u_equ,k) = x_s_u_fre;

                if isnan(sum(x_s_u_fre)),  disp('LS error'); end

            end  % N_k end

            x_nc_block(:,:,sblk) = x_nc_subblock;

            % demapping symbol in frequency domain  for each active user
            for u = 1:N_a

                data_u = (x_nc_subblock(u,:)).';
                mask_pos_sc_u = mask_pos_all(:,est_set(u));
                data_u_td = data_u(mask_pos_sc_u,:);  % time domain info
                data_u_fd =  Fo_DFT*data_u_td; % frequency domain info
                data_u_llr = qamdemod(data_u_fd, M_ord, 'gray','OutputType', 'llr', 'UnitAveragePower', true); % in bit dim
                data_recv_info((sblk-1)*N_info*N_bits+1:sblk*N_info*N_bits,blk,u,np) = data_u_llr;


            end

        end  % N_subblock end


    end  % N_block end
end % N_p end

%% channel decoding
% de-FFT, OFDM - finished during detection
% de-spread, s - finished during detection
% de-DFT, - finished during detection
% channel decoding
data_recv_info_all = reshape(data_recv_info,[N_info_pblock_bit*N_block,N_a,N_p]);

for u = 1:N_a
   
    u_idx = est_set(u);
    data_u_p_llr_all = data_recv_info_all(:,u,:);
    data_u_p_llr  =  log( mean(exp(data_u_p_llr_all),3) ); %  average decoding for distributed system
    % data_u_p_llr  =  data_u_p_llr_all(:,:,1); %  average decoding for distributed system
    switch Coded_Sel
        case 'LDPC'
            data_recv_info_bit(:,u_idx) = (ldpcDecode(data_u_p_llr', LDPC)).';
        case 'No'
            data_recv_info_bit(:,u_idx) = data_u_p_llr'<0;
    end
    data_recv_info_bit_Uncoded(:,u_idx) = data_u_p_llr'<0;
    % bit = double(data_u_p_llr<0);
    % errbits = sum(xor(bit,Coded_Info_bits(:,1)));
end


end  % function end

