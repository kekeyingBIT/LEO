function [rec_sig, Nvar, Agc] = InfoTrans(sframe,Hd_mtx,act_set,L_max,N_r, N_p,snr)
% transmit information 
N_a = length(act_set);

[data_len, N_s, N_u] = size(sframe);
N_BS = N_r(1)*N_r(2);
rec_sig = zeros(data_len+L_max-1,N_BS,N_p);

for u = 1:N_a
    
    u_idx = act_set(u);
    for s = 1:N_s
        h_u = Hd_mtx((u_idx-1)*L_max+1:u_idx*L_max,:,s);
        sig_u = sframe(:,s,u_idx);
        for n = 1:N_BS
            rec_sig(:,n,s) = rec_sig(:,n,s) + conv(sig_u,h_u(:,n));
        end  % Nr end

    end  % Ns end

end  % Na end

% add noise at differenct receiver separately & AGC
[nL,nR,nP] = size(rec_sig);
Nvar = zeros(1,nP);
Agc = zeros(1,nP);
for np = 1:N_p 
    % add noise
    sigpow = norm(rec_sig(:,:,np),'fro')^2/(nL*nR);
    nvar = 10^(-snr/10)*sigpow;
    rec_sig_np = rec_sig(:,:,np);
    rec_sig(:,:,np) = rec_sig(:,:,np) + ...
        sqrt(nvar/2)* ( randn(nL,nR) + 1i*randn(nL,nR) );
    % acquire agc
    tmp_sig  = [ real(rec_sig_np(:)), imag(rec_sig_np(:)) ];
    agc = 1/max(tmp_sig(:));
    
    % store para
    Nvar(np) = nvar;
    Agc(np) = agc;

end


end

