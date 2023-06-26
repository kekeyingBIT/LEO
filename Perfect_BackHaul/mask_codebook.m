function [mask_seq] = mask_codebook(sf)
% generate random mask sequence to decorrelation

switch  sf
    case -1 % random generated codebook
        sf = 768;
        cb_size = 100;
        nzn = sf/2;
        mask_seq.set = zeros(sf,cb_size);
        for cw = 1:cb_size
            codeword = zeros(sf,1);
            nzn_idx = randperm(sf,nzn);
            codeword(nzn_idx,:) = ones(length(nzn_idx),1);
            mask_seq.set(:,cw) =  codeword;
        end
        mask_seq.cr = nzn/sf;  
        
    otherwise  % support 1/2/3/4/5/6/9/10
        mask_seq.set = eye(sf);
        nzn = 1;
        mask_seq.cr = nzn/sf;
end
    
end

