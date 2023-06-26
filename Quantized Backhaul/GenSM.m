function [Phi] = GenSM(TS_symbol,n_IBI,n_ts,N_block,L_max,N_u)
% Generate sensing matrix 
Phi_u =  zeros(n_IBI,L_max);
Phi = zeros(n_IBI,L_max*N_u,N_block);

for sym = 1:N_block
    tr_ts = TS_symbol((sym-1)*n_ts+1:sym*n_ts,:);
    for uu = 1:N_u
        seq  = tr_ts(:,uu);
        for ll = 1:n_IBI
            win = seq(ll:ll + L_max-1);
            Phi_u(ll,:) = fliplr(win.');
        end
        Phi(:,(uu-1)*L_max +1 : uu*L_max,sym) = Phi_u;
    end
end

end

