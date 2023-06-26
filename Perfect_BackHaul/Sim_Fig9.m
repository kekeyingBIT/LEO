% code for rl
% investigate the impact of number of path number on the
% overall performance, same as sim3 but different in parameters
% ch.N_c = 2
%%
clc; clear; close all
rng(20230203)
diary Sim4Diary
warning off
addpath(genpath( './ldpcEnc'))
addpath(genpath( './csAlg'))
%% parameter setup
% simulation para
Rho = 0.08;      % observation ratio >=60
SNR = 12;                                    % typical SNR at satellite side
nSim = 1000;                                 % Simulation times

% User Terminal
user.N_u = 100;                              % total user number in served area
user.PaAll = [0.15:0.05:0.45];                             % active ratio among all the potential users
user.N_s = 3;                                % only valid for diversity transmissoin, initialize data stream number at the user side

% LEO config
leo.Height =550;                             % height in km
leo.Dis = 500;                               % distance  between two LEO
leo.N_r = [10,10];
leo.N_p = 3;                                 % initialize satellite number

% frequency/time domain parameters
fq.N_k = 540;                                % subcarrier number for each TSP frame
fq.scs = 15e3;                               % subcarrier interval
fq.Bw = fq.N_k*fq.scs;                       % system bandwidth

Tb = 1/fq.scs;                               % physical time duration for each TSP block
Ts = 1/fq.Bw;                                % physical time duration for each TSP symbol

% channel parameter
ch.f_c  = 14.5e9;                            % uplink carrier frequency
ch.N_c = 2;                                  % path number = 1, due the beamforming capability at tx
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

% memory allocation and scheme choice
CE_Sel = char('SOMP','Oracle-LS','OAMP-1','OAMP-2');         % choose JADCE algorithm
DD_Sel = char('Distributed','Centralized'); %                            % choose DD algorithm
CIR_Sel = char('SOMP','Oracle-LS','OAMP-1','OAMP-2');                                                  % choose CIR used in DD
DT = char('diversity');                                                    % transmission mode selection
Coded_Sel = char('No');  %   char('LDPC');                                              % Channel coding scheme, No/LDPC

Nmse = zeros(size(CE_Sel,1),length(user.PaAll ));
Acc  = zeros(size(CE_Sel,1),length(user.PaAll ));
Acc2  = zeros(size(CE_Sel,1),length(user.PaAll ));
Pe =  zeros(size(CE_Sel,1),length(user.PaAll ));
Pe2 =  zeros(size(CE_Sel,1),length(user.PaAll ));
Ber_coded  = zeros(size(DD_Sel,1),size(CIR_Sel,1),length(user.PaAll ));                 % [Symstem Model, Observations, CIR Info]
Ber_uncoded  = zeros(size(DD_Sel,1),size(CIR_Sel,1),length(user.PaAll ));                 % [Symstem Model, Observations, CIR Info]
%% Simulation Begin
t_all = 0; % timer
for sim = 1:nSim
    % pre-allocate memory
    HD_est = zeros(ch.L_max*user.N_u, prod(leo.N_r),leo.N_p,size(CE_Sel,1), length(user.PaAll));
    ACT_idx = zeros(user.N_u,size(CE_Sel,1),length(user.PaAll));
    Rec_Sig = cell(length(user.PaAll),1);

    for na = 1:length(user.PaAll)
        tic
        user.P_a = user.PaAll(na);
        % generate active user set
        act_set = sort(randperm(user.N_u, double(int32(user.N_u*user.P_a))),'ascend');  % active user set
        act_idx = zeros(user.N_u,1); act_idx(act_set) = ones(length(act_set),1);         % activity indicater

        % data generation for active users
        [InfoBitsCoded, InfoBitsSource,enc] = DataBitsStream(user,act_set,Coded_Sel,'diversity');

        % TS length & Length of nIBI, block number
        l_ts = L_ts; l_nIBI = L_nIBI;
        N_block = size(InfoBitsCoded,1)/log2(M_ord)/fq.N_k;
        % generata TS and sensing matrix
        TS_symbol = GenTS('iidg',l_ts,l_nIBI,N_block,user.N_u);
        TS_Sym  = TS_symbol;
        Phi = GenSM(TS_symbol,l_nIBI,l_ts,N_block,ch.L_max,user.N_u);
       
        % introduce mask sequence to decorrelation
        mask = mask_codebook(sf);
        % generate transmission frame
        [sframe,sframe_symbol,mask_pos_all,fr] = GenFrame(act_set,TS_symbol,InfoBitsCoded,N_block,M_ord,mask,fq);
        Fr = fr;


        % generate TSL channel -> assume UT side beamforming
        [Hd_mtx, Hf_mtx, Hd, Hf, ch_para] = TSL2(leo,user,ch,fq,act_set);
        
        % printf CE dimension
        fprintf('ob_ratio = %.4f, sp_ratio = %.4f, signal_dim = %d, sp_level = %d, ob_length = %d \n',...
            l_nIBI/size(Phi,2), ch.N_c*length(act_set)/size(Phi,2), size(Phi,2),ch.N_c*length(act_set),l_nIBI);

        % data frame transmision (already add noise)

        [rec_sig, Nvar, Agc] = InfoTrans(sframe,Hd_mtx,act_set,ch.L_max,leo.N_r,leo.N_p, SNR);
        % received pilot sequence
        ob_nIBI = ParsePilot(rec_sig,N_block,l_nIBI,leo,fr); % [nIBI, Nr, N_block, N_p]
        % save rec_sig for data detection
        Rec_Sig{na,1} = rec_sig;

        % JADCE
        for ce_sel = 1:size(CE_Sel,1)  % select JADCE algorithm
            ce_choice = deblank(CE_Sel(ce_sel,:));
            [Hd_est, est_idx,est_idx_init] = MultiLEO_JADCE(ce_choice,ob_nIBI,Phi,Nvar,Agc, ...
                Hd_mtx,act_idx,ch,leo,user,na,ch_para);
            Nmse(ce_sel,na) = Nmse(ce_sel,na) + ( norm(Hd_est(:)-Hd_mtx(:),2) / norm(Hd_mtx(:),2) )^2;
            Acc(ce_sel,na) = Acc(ce_sel,na) + leo.N_p*( sum( xor( act_idx,est_idx ) ) == 0 ); % Centralized
            Acc2(ce_sel,na) = Acc2(ce_sel,na) + sum( sum( xor( act_idx,est_idx_init ) ) == 0 ); % Distributed

            Pe(ce_sel,na) = Pe(ce_sel,na) + leo.N_p*sum( xor( act_idx,est_idx ) ); % Centralized
            Pe2(ce_sel,na) = Pe2(ce_sel,na) + sum( xor( act_idx,est_idx_init ),'all' ); % Distributed


            fprintf('JADCE: nSim = %d/%d, Alg = %s, #Ka = %d,  NMSE = %.4f dB, Acc = %.4f, Pe = %.4f, Pe2 = %.4f\n',...
                sim, nSim, ce_choice, user.PaAll(na)*user.N_u, ...
                10*log10(Nmse(ce_sel,na)/sim), ...
                Acc(ce_sel,na)/sim/leo.N_p, ...
                Pe(ce_sel,na)/user.N_u/sim/leo.N_p,...
                Pe2(ce_sel,na)/user.N_u/sim/leo.N_p);

            HD_est(:,:,:,ce_sel,na) = Hd_est;
            ACT_idx(:,ce_sel,na) = est_idx;
        end  % ce_sel end

        %     end  % na end
        %     % toc
        %
        %     % data detection
        %     for na = 1:length(user.PaAll)
        TS_symbol = TS_Sym;
        fr = Fr;

        for cir_sel = 1:size(CIR_Sel,1)  % select CIR algorithm
            cir_choice = deblank(CIR_Sel(cir_sel,:));

            % mapping from cir_sel-> ce_sel
            cir_str = CIR_Sel(cir_sel,:);
            for cp = 1:size(CE_Sel,1)
                if strcmp(deblank(cir_str),deblank(CE_Sel(cp,:)))
                    ce_idx_sel = cp;
                end
            end

            Hd_est_tmp = HD_est(:,:,:,ce_idx_sel,na);
            act_est_tmp = ACT_idx(:,ce_idx_sel,na);
            rec_sig = Rec_Sig{na,1};
            % prepossing
            [PilotIBI] = calPilotIBI(TS_symbol,Hd_est_tmp,act_est_tmp,N_block,fr);
            [chMtx1] = GenFreqCh(Hd_est_tmp,act_est_tmp,mask_pos_all,fr);

            para.LDPC = enc; para.mask = mask;  para.fr = fr;
            para.N_block = N_block; para.L_max = ch.L_max; para.M_ord = M_ord; para.Nvar = Nvar;

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
                end

                % save result
                switch DT
                    case 'diversity'
                        data_trans_info_bit = squeeze( InfoBitsSource(:,1,:) );
                        data_trans_info_bit_raw = squeeze( InfoBitsCoded(:,1,:) );

                        NumErrBitsSource = sum(sum(data_recv_info_bit(:,act_set) ~= data_trans_info_bit(:,act_set) ) );
                        NumErrBitsRaw = sum(sum(data_recv_info_bit_Uncoded(:,act_set) ~= data_trans_info_bit_raw(:,act_set) ) );


                        Ber_coded(dd_sel,cir_sel,na)  =  Ber_coded(dd_sel,cir_sel,na)  +  NumErrBitsSource;  % Coded Ber
                        Ber_uncoded(dd_sel,cir_sel,na)  =  Ber_uncoded(dd_sel,cir_sel,na)  +  NumErrBitsRaw;  % Uncoded Ber



                    case 'multiplex'
                        error('Not supported for multiplex transmission')
                end

                % print status information
                fprintf('DD: nSim = %d/%d, CIR = %s, DD = %s, #Ka = %d, Coded Ber = %.6f, Uncoded Ber = %.6f\n',...
                    sim, nSim, cir_choice, dd_choice, user.PaAll(na)*user.N_u, ...
                    Ber_coded(dd_sel,cir_sel,na)/numel(InfoBitsSource(:,act_set))/sim,...
                    Ber_uncoded(dd_sel,cir_sel,na)/numel(InfoBitsCoded(:,act_set))/sim );
            end  % dd_sel end
        end % ce_sel end

    end % na end

    % save simulation result
    BerCoded_th = Ber_coded/numel(InfoBitsSource(:,act_set))/sim;
    BerUncoded_th = Ber_uncoded/numel(InfoBitsCoded(:,act_set))/sim;
    Nmse_th = 10*log10(Nmse/sim);
    Acc_th =  Acc/sim*100/leo.N_p;  % convert to 100%
    Acc_th2 =  Acc2/sim*100/leo.N_p;  % convert to 100%

    Pe_th = Pe/user.N_u/sim/leo.N_p;
    Pe_th2 = Pe2/user.N_u/sim/leo.N_p;

    data.BerCoded_th = BerCoded_th;
    data.BerUncoded_th = BerUncoded_th;
    data.Nmse_th = Nmse_th;
    data.Acc_th = Acc_th;
    data.Acc_th2 = Acc_th2;

    data.Pe_th = Pe_th;
    data.Pe_th2 = Pe_th2;

    % save raw data for recovery
    data_raw.BerCoded_th =  Ber_coded;
    data_raw.BerUncoded_th =  Ber_uncoded;
    data_raw.Nmse = Nmse;
    data_raw.Acc = Acc;
    data_raw.Acc2 = Acc2;
    data_raw.Pe = Pe;
    data_raw.Pe2 = Pe2;

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

