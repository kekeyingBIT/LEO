%% version information
% Fig 13, add quantization to feedback
% BER_vs_pilot length, LMMSE for Alg. 3 rather than LS
%%
clc; clear; close all
rng(20221019)
diary mydiary
warning off
addpath(genpath( './ldpcEnc'))
addpath(genpath( './csAlg'))
%% parameter setup
% simulation para
Rho = [0.036,0.04:0.01:0.12,0.15];      % observation ratio >=60
SNR = 12;%[5,8,12,15];                                    % typical SNR at satellite side
nSim = 500;                                 % Simulation times

% User Terminal
user.N_u = 100;                              % total user number in served area
user.P_a = 0.15;                             % active ratio among all the potential users
user.N_s = 3;                                % data stream number at the user side

% LEO config
leo.Height =550;                             % height in km
leo.Dis = 500;                               % distance  between two LEO
leo.N_r = [10,10]; % [5,5]                           % receive antenna size
leo.N_p = 3;                                 % satellite number

% frequency/time domain parameters
fq.N_k = 540;                                % subcarrier number for each TSP frame
fq.scs = 15e3;                               % subcarrier interval
fq.Bw = fq.N_k*fq.scs;                       % system bandwidth

Tb = 1/fq.scs;                               % physical time duration for each TSP block
Ts = 1/fq.Bw;                                % physical time duration for each TSP symbol

% channel parameter
ch.f_c  = 14.5e9;                            % uplink carrier frequency
ch.N_c = 3;                                  % path number = 1, due the beamforming capability at tx
ch.L_max = 17;                               % tap number of maximum diffential delay 
ch.tau_max = ch.L_max*Ts;

ch.plos = 1;                               % Los probability of each link for urban scenario
ch.Kf = 10;                                  % Rician factor in dB 
ch.ag = 0;                                  % array gain of beamforming in dB

% length parameter
L_nIBI = double(int32(ch.L_max*user.N_u*Rho));   % length of no IBI region
L_ts = double(L_nIBI+ch.L_max-1);            % length of TS
L_all = L_ts + fq.N_k;                          % total length of a TSP frame

T_nIBI = double(L_nIBI)*Ts;                  % time duration of no IBI region
T_ts = double(L_ts)*Ts;                      % time duration of TS
T_all = L_all*Ts;                            % time duration a TSP frame

% data transmission parameter
M_ord = 4;                                   % modulation order
N_bits = log2(M_ord);                        % information bits of per symbol
% e_norm = (M_ord==2)+...
%     (M_ord~=2)*sqrt((M_ord-1)/6*(2^2));      % energy normalized factor

sf = 1;                                      % lenght of spreading factor
% quantization feed bits
Qbits = 3;

% memory allocation and scheme choice
CE_Sel =  char('OAMP-2');         % choose JADCE algorithm
DD_Sel = char('Centralized_Quant_LS', 'Centralized_Quant_LS_mix','Centralized_Quant_VI_2','Centralized_Quant_VI_1','Centralized_Unquant_LS');                                % choose DD algorithm
% DD_Sel = char('Centralized_Quant_LS','Centralized_Quant_VI_3','Centralized_Quant_VI_2','Centralized_Quant_VI_1','Centralized_Quant_VI_0','Centralized_Unquant_LS');                                % choose DD algorithm
% DD_Sel = char('Centralized_Quant_LS','Centralized_Quant_VI_2','Centralized_Unquant_LS');                                % choose DD algorithm
CIR_Sel = char('OAMP-2');                                  % choose CIR used in DD
DT = char('diversity');                                                    % transmission mode selection
Coded_Sel = char('No');                                                  % Channel coding scheme, No/LDPC

Nmse = zeros(size(CE_Sel,1),length(Rho),length(SNR));
Acc  = zeros(size(CE_Sel,1),length(Rho),length(SNR));
Pe =  zeros(size(CE_Sel,1),length(Rho),length(SNR));
Ber_coded  = zeros(size(DD_Sel,1),length(Rho),size(CIR_Sel,1),length(SNR));                 % [Symstem Model, Observations, CIR Info]
Ber_uncoded  = zeros(size(DD_Sel,1),length(Rho),size(CIR_Sel,1),length(SNR));                 % [Symstem Model, Observations, CIR Info]
%% Simulation Begin
t_all = 0; % timer
for sim = 1:nSim
    tic
    % generate active user set
    act_set = sort(randperm(user.N_u, double(int32(user.N_u*user.P_a))),'ascend');  % active user set 
    act_idx = zeros(user.N_u,1); act_idx(act_set) = ones(length(act_set),1);         % activity indicater
    
    % generate TSL channel -> assume UT side beamforming
    [Hd_mtx, Hf_mtx, Hd, Hf, ch_para] = TSL(leo,user,ch,fq,act_set);
        
    % data generation for active users
    [InfoBitsCoded, InfoBitsSource,enc] = DataBitsStream(user,act_set,Coded_Sel,'diversity');
    
    % pre-allocate memory 
    HD_est = zeros(ch.L_max*user.N_u, prod(leo.N_r),leo.N_p,size(CE_Sel,1), length(Rho),length(SNR));
    ACT_idx = zeros(user.N_u,size(CE_Sel,1),length(Rho),length(SNR));
    
    % fr_len =  size(TS_symbol,1) + size(InfoBitsCoded,1)/N_bits/mask.cr;
    % Rec_Sig = zeros(fr_len,numel(leo.N_r),leo.N_p,length(SNR));
    Rec_Sig = cell(length(Rho),length(SNR));
    TS_Sym = cell(length(Rho),1);
    Fr = cell(length(Rho),1);
    
    for ob = 1:length(Rho)

        % TS length & Length of nIBI, block number
        l_ts = L_ts(ob); l_nIBI = L_nIBI(ob);
        N_block = size(InfoBitsCoded,1)/log2(M_ord)/fq.N_k; 
        % generata TS and sensing matrix
        TS_symbol = GenTS('iidg',l_ts,l_nIBI,N_block,user.N_u); 
        TS_Sym{ob}  = TS_symbol;
        Phi = GenSM(TS_symbol,l_nIBI,l_ts,N_block,ch.L_max,user.N_u);
        % printf CE dimension
        fprintf('ob_ratio = %.4f, sp_ratio = %.4f, signal_dim = %d, sp_level = %d, ob_length = %d \n',...
            l_nIBI/size(Phi,2), ch.N_c*length(act_set)/size(Phi,2), size(Phi,2),ch.N_c*length(act_set),l_nIBI);
        
        % introduce mask sequence to decorrelation
        mask = mask_codebook(sf);
        % generate transmission frame
        [sframe,sframe_symbol,mask_pos_all,fr] = GenFrame(act_set,TS_symbol,InfoBitsCoded,N_block,M_ord,mask,fq);
        Fr{ob} = fr;
        
        NVAR = zeros(length(SNR),leo.N_p);
        for snr = 1:length(SNR)
            % data frame transmision (already add noise)
            [rec_sig, Nvar, Agc] = InfoTrans(sframe,Hd_mtx,act_set,ch.L_max,leo.N_r,leo.N_p, SNR(snr));
            NVAR(snr,:) = Nvar;
            % received pilot sequence
            ob_nIBI = ParsePilot(rec_sig,N_block,l_nIBI,leo,fr); % [nIBI, Nr, N_block, N_p]
            % save rec_sig for data detection
            Rec_Sig{ob,snr} = rec_sig;
            
            % JADCE
            for ce_sel = 1:size(CE_Sel,1)  % select JADCE algorithm
                ce_choice = deblank(CE_Sel(ce_sel,:));
                [Hd_est, est_idx] = MultiLEO_JADCE_SNR(ce_choice,ob_nIBI,Phi,Nvar,Agc, ...
                    Hd_mtx,act_idx,ch,leo,user,ob,snr,ch_para);
                Nmse(ce_sel,ob,snr) = Nmse(ce_sel,ob,snr) + ( norm(Hd_est(:)-Hd_mtx(:),2) / norm(Hd_mtx(:),2) )^2;
                Acc(ce_sel,ob,snr) = Acc(ce_sel,ob,snr) + ( sum( xor( act_idx,est_idx ) ) == 0 );
                Pe(ce_sel,ob,snr) = Pe(ce_sel,ob,snr) + sum( xor( act_idx,est_idx ) );
                
                fprintf('JADCE: nSim = %d/%d, Alg = %s, rho = %.3f, snr = %d dB, NMSE = %.4f dB, Acc = %.4f, Pe = %.4f\n',...
                    sim, nSim, ce_choice, Rho(ob), SNR(snr),...
                    10*log10(Nmse(ce_sel,ob,snr)/sim), Acc(ce_sel,ob,snr)/sim, Pe(ce_sel,ob,snr)/user.N_u/sim);
                
                HD_est (:,:,:,ce_sel,ob,snr) = Hd_est;
                ACT_idx(:,ce_sel,ob,snr) = est_idx;  
            end  % ce_sel end
            
        end   % snr end
    end  % ob end
    % toc

    % data detection
    for ob = 1:length(Rho)
        TS_symbol = TS_Sym{ob};
        fr = Fr{ob};
        
        for cir_sel = 1:size(CIR_Sel,1)  % select CIR algorithm
            cir_choice = deblank(CIR_Sel(cir_sel,:));
            
            % mapping from cir_sel-> ce_sel
            cir_str = CIR_Sel(cir_sel,:);
            for cp = 1:size(CE_Sel,1) 
                if strcmp(deblank(cir_str),deblank(CE_Sel(cp,:)))
                    ce_idx_sel = cp;
                end
            end
            % ce_idx_sel = 2;
            
            for snr = 1:length(SNR)
                Hd_est_tmp = HD_est(:,:,:,ce_idx_sel,ob,snr);
                act_est_tmp =ACT_idx(:,ce_idx_sel,ob,snr);
                rec_sig = Rec_Sig{ob,snr};
                % prepossing
                [PilotIBI] = calPilotIBI(TS_symbol,Hd_est_tmp,act_est_tmp,N_block,fr);
                [chMtx1] = GenFreqCh(Hd_est_tmp,act_est_tmp,mask_pos_all,fr);
                
                para.LDPC = enc; para.mask = mask;  para.fr = fr;
                para.N_block = N_block; para.L_max = ch.L_max; para.M_ord = M_ord; 
                para.Nvar = NVAR(snr,:);
                
                for dd_sel = 1:size(DD_Sel,1) % select centralized or distributed processing

                    dd_choice = deblank(DD_Sel(dd_sel,:));
                    
                    % data detection
                    switch dd_choice
                        case 'Distributed'
                            [data_recv_info_bit,data_recv_info_bit_Uncoded] = DD(rec_sig,PilotIBI,chMtx1,act_est_tmp,mask_pos_all,dd_choice,para,Coded_Sel,...
                                InfoBitsCoded,InfoBitsSource,sframe_symbol);
                        case 'Centralized'
                            [data_recv_info_bit,data_recv_info_bit_Uncoded] = DD2(rec_sig,PilotIBI,chMtx1,act_est_tmp,mask_pos_all,dd_choice,para,Coded_Sel,...
                                InfoBitsCoded,InfoBitsSource,sframe_symbol);
                        case { 'Centralized_Quant_LS_mix','Centralized_Quant_VI_3','Centralized_Quant_VI_2','Centralized_Quant_VI_1','Centralized_Quant_VI_0','Centralized_Quant_LS','Centralized_Unquant_LS'}  % centralized detection with quantized feedback
                            [data_recv_info_bit,data_recv_info_bit_Uncoded] = DD3(rec_sig,PilotIBI,chMtx1,act_est_tmp,mask_pos_all,dd_choice,para,Coded_Sel,...
                                InfoBitsCoded,InfoBitsSource,sframe_symbol,Qbits);
                            
                    end
                    
                    % save result
                    switch DT
                        case 'diversity'
                            data_trans_info_bit = squeeze( InfoBitsSource(:,1,:) );
                            data_trans_info_bit_raw = squeeze( InfoBitsCoded(:,1,:) );
                            
                            NumErrBitsSource = sum(sum(data_recv_info_bit(:,act_set) ~= data_trans_info_bit(:,act_set) ) );
                            NumErrBitsRaw = sum(sum(data_recv_info_bit_Uncoded(:,act_set) ~= data_trans_info_bit_raw(:,act_set) ) );
                            
                            
                            Ber_coded(dd_sel,ob,cir_sel,snr)  =  Ber_coded(dd_sel,ob,cir_sel,snr)  +  NumErrBitsSource;  % Coded Ber
                            Ber_uncoded(dd_sel,ob,cir_sel,snr)  =  Ber_uncoded(dd_sel,ob,cir_sel,snr)  +  NumErrBitsRaw;  % Uncoded Ber
                            
                            
                            
                        case 'multiplex'
                            error('Not supported for multiplex transmission')
                    end
                    
                    % print status information
                    fprintf('DD: nSim = %d/%d, CIR = %s, DD = %s, rho = %.3f, snr = %d, Coded Ber = %.6f, Uncoded Ber = %.6f\n',...
                        sim, nSim, cir_choice, dd_choice, Rho(ob), SNR(snr),...
                        Ber_coded(dd_sel,ob,cir_sel,snr)/numel(InfoBitsSource(:,act_set))/sim,...
                        Ber_uncoded(dd_sel,ob,cir_sel,snr)/numel(InfoBitsCoded(:,act_set))/sim );
                end  % dd_sel end

            end % snr end
        end % ce_sel end

    end % ob end

   % save simulation result
    BerCoded_th = Ber_coded/numel(InfoBitsSource(:,act_set))/sim;
    BerUncoded_th = Ber_uncoded/numel(InfoBitsCoded(:,act_set))/sim;
    Nmse_th = 10*log10(Nmse/sim);
    Acc_th =  Acc/sim*100;  % convert to 100%
    Pe_th = Pe/user.N_u/sim;

    data.BerCoded_th = BerCoded_th; 
    data.BerUncoded_th = BerUncoded_th;
    data.Nmse_th = Nmse_th;
    data.Acc_th = Acc_th;
    data.Pe_th = Pe_th;
    
    % save raw data for recovery
    data_raw.BerCoded_th =  Ber_coded;
    data_raw.BerUncoded_th =  Ber_uncoded;
    data_raw.Nmse = Nmse;
    data_raw.Acc = Acc;
    data_raw.Pe = Pe;
    
    % save temporal data 
    if mod(sim,1) == 0
        data_name = ['./Data/',num2str(sim),'th','.mat'];
        save(data_name,'data');
        
        data_raw_name = ['./Data/','raw',num2str(sim),'th','.mat'];
        save(data_raw_name,'data_raw');
        
    end
    
    
    
    t_end = toc;
    t_all = t_all + t_end;
    fprintf('Time has passed by: %.2f s @ %d -th sim. \n', t_all, sim);

end  % sim end

%% plot simulation result

