function [TS_symbol] = GenTS(ts_sel,n_ts,n_nIBI,N_block,N_u)
% training sequence selection, all data stream/LEO share the same TS
switch ts_sel
    case 'iidg'
         TS_symbol = sqrt(1/2)*(randn(n_ts*N_block,N_u) + 1i*randn(n_ts*N_block,N_u));
    case 'binary'
         TS_symbol = 2*random('Binomial',1,0.5,[n_ts*N_block,N_u])-1;
%     case 'binary'
%          TS_symbol = zeros(n_ts*N_block, N_u);       
%          % pn sequence: m seq
%          for u = 1:N_u
%              h1 = commsrc.pn('GenPoly',     [1 1 0 1 1 0 0 0 0 0 0 0 0 1],...%[13 12 10 9 0] ;
%                    'InitialStates',[0 0 0 0 0 0 0 0 0 0 0 0 1],...  
%                    'CurrentStates',[0 0 0 0 0 0 0 0 0 0 0 0 1],...  ?
%                    'shift',0,...                                   
%                    'NumBitsOut',n_ts*N_block);                            
%             PN_seq = -generate(h1).*2+1; 
%             TS_symbol(:,u) = PN_seq;
%          end
    case 'qam'    
        % iid symbol drawn from qam constellation set
        DACReso = 2*[1 2 3]; % 1/2/3 bit DAC
        bit= DACReso(1); % 4QAM 导频/16QAM/64QAM导频
        Morder = 2^bit; % 调制阶数 4QAM/16QAM/64QAM
        %             pilot_bit = randi([0 1],1,bit*L*N); % bit generation
        pilot_bit = randi([0 1],1,bit*n_ts*N_block*N_u); % bit generation
        bitmatrix = reshape(pilot_bit,length(pilot_bit)/bit,bit);
        symbols = bi2de(bitmatrix);
        symbol = qammod(symbols,Morder,'gray').';
        Scale = modnorm(symbol,'avpow',1); % normalization
        seq_all = symbol*Scale;
        TS_symbol= reshape(seq_all,n_ts*N_block,N_u);  % N 用户每符号功率和为1，总为L个导频能量

end

