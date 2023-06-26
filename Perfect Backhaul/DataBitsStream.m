function [InfoBitsCoded, InfoBitsUnCoded,enc] = DataBitsStream(user,act_set,Enc_Sel,DT)

switch Enc_Sel

    case 'LDPC'
        % ldpc parameter
        blkSize = 1024; codeRate = '1/2';  % which generate 2112 uncoded bits and 4320 coded bits
        LDPC = ldpcGet(blkSize, codeRate);
        enc = LDPC;
        % data stream
        InfoBitsCoded = zeros(LDPC.numTotBits,user.N_s, user.N_u);
        InfoBitsUnCoded = zeros(LDPC.numInfBits,user.N_s,user.N_u);
        % choose tranmission mode
        switch DT
            % different rf chain share the same data 
            case 'diversity' 
                for u = 1:length(act_set)
                    data = randi([0 1], 1, LDPC.numInfBits);
                    dataEnc = ldpcEncode(data, LDPC);
                    InfoBitsUnCoded(:,:,act_set(u)) = repmat(data.',[1,user.N_s,1]);
                    InfoBitsCoded(:,:,act_set(u)) = repmat(dataEnc.',[1,user.N_s,1]);
                end
            % different rf chain transmit different data 
            case 'multiplex'
                for u = 1:length(act_set)
                    for s = 1:user.N_s
                        data = randi([0 1], 1, LDPC.numInfBits);
                        dataEnc = ldpcEncode(data, LDPC);
                        InfoBitsUnCoded(:,s,act_set(u)) = data.';
                        InfoBitsCoded(:,s,act_set(u)) = dataEnc.';
                    end
                end
            otherwise
                error ('Unsupported data transmission method');
        end 

    case 'No'
        % ldpc parameter
        blkSize = 1024; codeRate = '1/2';  % which generate 2112 uncoded bits and 4320 coded bits
        LDPC = ldpcGet(blkSize, codeRate);
        enc = LDPC;
        % data stream
        InfoBitsCoded = zeros(LDPC.numTotBits,user.N_s, user.N_u);
        InfoBitsUnCoded = zeros(LDPC.numTotBits,user.N_s,user.N_u);
        % choose tranmission mode
        switch DT
            % different rf chain share the same data
            case 'diversity'
                for u = 1:length(act_set)
                    data = randi([0 1], 1, LDPC.numTotBits);
                    InfoBitsUnCoded(:,:,act_set(u)) = repmat(data.',[1,user.N_s,1]);
                    InfoBitsCoded(:,:,act_set(u)) = repmat(data.',[1,user.N_s,1]);
                end
                % different rf chain transmit different data
            case 'multiplex'
                for u = 1:length(act_set)
                    for s = 1:user.N_s
                        data = randi([0 1], 1, LDPC.numTotBits);
                        InfoBitsUnCoded(:,s,act_set(u)) = data.';
                        InfoBitsCoded(:,s,act_set(u)) = data.';
                    end
                end
            otherwise
                error ('Unsupported data transmission method');
        end

    case 'Conv'
        error('Convolution coding is not supported now ');
    otherwise 
        error('Unsupported channel coding')

end


end

