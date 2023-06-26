function [Hd_est_refined,es] = chRebuild(Hd_est,DoA_est,est_idx,supp,ch,leo,user,ch_para, H_true)
%% parameters
N_r = leo.N_r;  Nx= N_r(1); Ny= N_r(2);
N_u = user.N_u;
L_max = ch.L_max; 
est_set = find(est_idx == 1);
%% active user DoA 
N_a = length(est_set);
Asv = zeros(Nx*Ny,N_a);
n_x = (0:Nx-1).';
n_y = (0:Ny-1).';
PHIs  = DoA_est(:,1).';
THETAs = DoA_est(:,2).';
e_x = exp(-1i*pi*n_x*  ( sind(PHIs).*cosd(THETAs) )  );
e_y = exp(-1i*pi*n_y*  ( sind(PHIs).*sind(THETAs) )  );
for u = 1: N_a
    Asv(:,u) = kron(e_y(:,u),e_x(:,u));
end
%% active user delay
tap = zeros(L_max,N_a);
for u = 1:N_a
    idx = est_set(u);
    tmp  = find ( supp( (idx-1)*L_max+1:idx*L_max) > 0 ) ;
    tap_u  = zeros(L_max,1);
    tap_u (tmp) = ones(length(tmp),1);
    tap(:,u) = tap_u;
end
%% active user gain
gain  = zeros(N_a,L_max);
for u =  1: N_a
    idx = est_set(u);
    seq_idx = (idx-1)*L_max+1:idx*L_max;
    supp_u  = find ( supp(seq_idx) > 0 )   ;
    Hk = Hd_est(seq_idx,: );
    hu_est = (Hk(supp_u,:));
    asv_u  = Asv(:,u);
    for tp = 1:length(supp_u)
        h_tap = hu_est(tp,:).';  % h_tap = asv_u*a
        gain(u,supp_u(tp)) = (asv_u'*asv_u)\asv_u'*h_tap;
    end
end


%% delay_channel_init

% Hd_est_refined = zeros(N_u*L_max,Nx*Ny);
Hd_est_refined = Hd_est; 

for u = 1:N_a
    u_idx = est_set(u);
    coef_u = gain(u,:).';
    ch_u = coef_u*Asv(:,u).';     
    seq_idx = (u_idx-1)*L_max+1:u_idx*L_max;

    % Hd_est_refined(seq_idx,:) = ch_u;

    supp_u = find(supp(seq_idx)>0);
    Hd_est_refined(seq_idx(supp_u),:) = ch_u(supp_u,:);

end 

es.phi = DoA_est(:,1);
es.theta = DoA_est(:,2);
es.ga = gain;
es.tap = tap;







end

