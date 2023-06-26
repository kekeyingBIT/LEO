function [Hd_mtx, Hf_mtx, Hd, Hf, ch_para] = TSL(leo,user,ch,fq,act_set)

% broadband MIMO channel for the delay domain

%% funvtion parameters
% Input parameters
% leo parameters:
    % leo.Height: LEO height
    % leo.Dis: distance  between two LEO
    % leo.N_r: LEO antenna number
    % leo.N_p: LEO satellte number

% user parameters:
    % user.N_u: potential user number
    % user.P_a: active probability
    % user.N_s: stream number for each user
% channel parameters
    % ch.f_c: carrier frequency
    % ch.N_c: multi-path number
    % ch.L_max: maximum  differential delay
    % ch.K: rician factor in dB
% frequency domain parameters
    % fq.N_k: subcarrier number
    % fq.scs: subcarier interval in Hz
    % fq.Bw: system bandwith
% act_set parameters
    % act_set: active user set

% output parameters
    % Hd_mtx: 2D  delay domain channel 
    % Hf_mtx: 2D frequency domain channel
    % Hd: 3D delay domain channel
    % Hf: 3D delay domain channel
    % ch_para: key parameters such as AoA & gain

%% generate users distribution and satellite network layout

% generate user distribution
N_a = length(act_set);

cellCenter = [0,0];
theta = linspace(0,2*pi,4)+pi/2;
radius = leo.Dis/sqrt(3);
diameter = leo.Dis/2;
hexaPointX = cellCenter(1) + radius.*cos(theta);
hexaPointY = cellCenter(2) + radius.*sin(theta);

refPoint = [0,0,leo.Height].';
numMU = 0;
TxPosition = zeros(3,N_a);
while numMU < N_a

    p = 2*diameter*rand(1,2)-diameter;
    flagIn = ( p(2)>hexaPointY(2) ) &&  ...
        ( p(2)< -sqrt(3)*p(1)+hexaPointY(1) ) && ...
        ( p(2)< sqrt(3)*p(1)+hexaPointY(1) );

    if numMU >= 1 
        pNow = [p,0].';
        LocExist = TxPosition(:,1:numMU);

        pNowVec = pNow-refPoint;
        pExistVec = LocExist - refPoint;

        pNowNorm = sqrt( sum(abs(pNowVec).^2,1) );
        pExistNorm = sqrt( sum(abs(pExistVec).^2,1) );

        CorrMat = pNowVec'*pExistVec;
        CorrNorm = pNowNorm*pExistNorm;
        Corr = CorrMat./CorrNorm;
        
        th_corr = 0.999;
        flagCorr = ( max(Corr(:)) > th_corr  );
    else
        flagCorr = 0;   % for the first user
    end

    if flagIn && ~flagCorr
        numMU = numMU +1;
        TxPosition(1:2,numMU) = p.';
    end       

end

% generate satellite distribution
RxPosition = zeros(3,leo.N_p);
RxPosition(1:2,1) =  [hexaPointX(1),hexaPointY(1)].';
RxPosition(1:2,2) =  [hexaPointX(2),hexaPointY(2)].';
RxPosition(1:2,3) =  [hexaPointX(3),hexaPointY(3)].';
RxPosition(3,:) = repmat(leo.Height,[1,leo.N_p]);

% visiualize UT and LEO distribution
% figure, grid on
% plot3(hexaPointX,hexaPointY,zeros(1,4),'b-', 'LineWidth',1.2); hold on
% scatter3(TxPosition(1,:), TxPosition(2,:),TxPosition(3,:),'ro', 'LineWidth',1.2);
% scatter3(RxPosition(1,:), RxPosition(2,:),RxPosition(3,:),'ko', 'LineWidth',1.2);

%% generate channel parameters

% steering vector
phis1 = 180/pi*atan( abs((TxPosition(3,:)-RxPosition(3,1))) ./ sqrt( (TxPosition(1,:)-RxPosition(1,1)).^2 + (TxPosition(2,:)-RxPosition(2,1)).^2 ) );
phis2 = 180/pi*atan( abs((TxPosition(3,:)-RxPosition(3,1))) ./ sqrt ( (TxPosition(1,:)-RxPosition(1,2)).^2 + (TxPosition(2,:)-RxPosition(2,2)).^2 ) );
phis3 = 180/pi*atan( abs((TxPosition(3,:)-RxPosition(3,1))) ./ sqrt( (TxPosition(1,:)-RxPosition(1,3)).^2 + (TxPosition(2,:)-RxPosition(2,3)).^2 ) );

thetas1 = 180/pi*atan2( (TxPosition(2,:)-RxPosition(2,1)), (TxPosition(1,:)-RxPosition(1,1)) );
thetas2 = 180/pi*atan2( (TxPosition(2,:)-RxPosition(2,2)), (TxPosition(1,:)-RxPosition(1,2)) );
thetas3 = 180/pi*atan2( (TxPosition(2,:)-RxPosition(2,3)), (TxPosition(1,:)-RxPosition(1,3)) );

ds1 = sqrt( (TxPosition(1,:)-RxPosition(1,1)).^2 + (TxPosition(2,:)-RxPosition(2,1)).^2 + (TxPosition(3,:)-RxPosition(3,1)).^2 );
ds2 = sqrt( (TxPosition(1,:)-RxPosition(1,2)).^2 + (TxPosition(2,:)-RxPosition(2,2)).^2 + (TxPosition(3,:)-RxPosition(3,2)).^2 );
ds3 = sqrt( (TxPosition(1,:)-RxPosition(1,3)).^2 + (TxPosition(2,:)-RxPosition(2,3)).^2 + (TxPosition(3,:)-RxPosition(3,3)).^2 );

PHIs = [phis1; phis2; phis3];
THETAs = [thetas1; thetas2; thetas3];
Ds = [ds1; ds2; ds3];

Asv = zeros(prod(leo.N_r),N_a,leo.N_p);
n_x= (0:leo.N_r(1)-1).';
n_y = (0:leo.N_r(2)-1).';
for np = 1:leo.N_p
    e_x = exp( -1i*pi*n_x* sind(PHIs(np,:)).*cosd(THETAs(np,:))   ); % [nx, nact];
    e_y = exp( -1i*pi*n_y* sind(PHIs(np,:)).*sind(THETAs(np,:))   ); % [nx, nact];
    for u = 1:N_a
        Asv(:,u,np) = kron(e_y(:,u),e_x(:,u));
    end
end

% equivalent channel gain after considering power control, path loss, analog beamforming, loss prob at tx
kf = 10^(ch.Kf/10);  % rician factor (original from channel self)
ag = 10^(ch.ag/10);  % array gain (original from analog beamforming)
los_idx = binornd(1,ch.plos,[N_a,user.N_s]); % assign los probability for each stream
init_gain = sqrt(1/2)*(randn(N_a,leo.N_p)+1i*randn(N_a,leo.N_p)); % randomness come from impefect power control and phase noise (equals 1 when perfect)
gain = zeros(N_a, ch.N_c, leo.N_p);
for u = 1: N_a
    for np = 1:leo.N_p
        if los_idx(u,np) == 1  % p-th stream channel of u-th UT is los
            los_gain = sqrt(kf/(kf+1))*init_gain(u,np); 
            nlos_gain = sqrt(1/(kf+1))*init_gain(u,np).* ( sqrt(1/2)*(randn(1,ch.N_c-1)+1i*randn(1,ch.N_c-1)) ); % randomness come from nlos
            gain(u,:,np) = [los_gain, nlos_gain];  % overall channel gain for all paths, [Na, Nc,Ns]
        elseif los_idx(u,np) == 0  % no los path between u-th UT and p-th LEO
            nlos_gain_sel  = sqrt(1/(kf+1))*init_gain(u,np).* ( sqrt(1/2)*(randn+1i*randn) );  % choose arbitrary NLOS for transmission
            nlos_gain_other  = sqrt(1/(kf+1))*init_gain(u,np).* ( sqrt(1/2)*(randn(1,ch.N_c-1)+1i*randn(1,ch.N_c-1)) );
            gain(u,:,np) = [nlos_gain_sel, nlos_gain_other];
        end % los prob end
    end % N_p end
end  % N_a end

% equivalent delay after considering timing advance
tap = zeros(ch.L_max, N_a, leo.N_p);
for np = 1:leo.N_p
    % randperm tap
    tap_tmp = zeros(ch.N_c,N_a);
    for u = 1:N_a
        tap_sel =randperm(ch.L_max,ch.N_c);
        tap_tmp(:, u) = tap_sel.';
    end   

    % alignment
    if ~ismember(1,tap_tmp)
        tap_tmp( randi(ch.N_c), randi(N_a) ) = 1;
    end
    if ~ismember(ch.L_max,tap_tmp)
        tap_tmp( randi(ch.N_c), randi(N_a) ) = ch.L_max;
    end
    tap_tmp = sort(tap_tmp,1,'ascend');
    
    % store tap
    for u = 1:N_a
        tap(tap_tmp(:,u),u,np) = ones(ch.N_c,1);
    end
end

%% generate channel 
% delay domain channel
Hd_mtx = zeros(user.N_u*ch.L_max,prod(leo.N_r),leo.N_p); % [Nu*L, Nr, Np]
Hd = zeros(ch.L_max,prod(leo.N_r),user.N_u,leo.N_p);     % [L, Nr, Nu, Np]
for np = 1:leo.N_p
    for u = 1:N_a
        
        tap_u = tap(:,u,np);                        % delay
        coef_u = gain(u,:,np).';                       % gain
        ch_u = zeros(ch.L_max,prod(leo.N_r));      % [L, Nr]
        ch_u( tap_u==1,:) = coef_u * Asv(:,u,np).';
        
        u_idx = act_set (u);                        % user id
        Hd_mtx((u_idx-1)*ch.L_max+1:u_idx*ch.L_max,:,np) = ch_u;
        Hd(:,:,u_idx,np) = ch_u;
    end
end

% frequency domain channel
Hf_mtx = zeros(user.N_u*fq.N_k, prod(leo.N_r), leo.N_p );  % [Nu*Nk, Nr, Np]
Hf = zeros(fq.N_k, prod(leo.N_r),user.N_u, leo.N_p);    % [Nk, Nr, Nu, Np]

for np = 1:leo.N_p
    for u = 1:N_a
        tap_u = tap(:,u,np); 
        tap_supp = find(tap_u == 1);
        ch_u = zeros(fq.N_k, prod(leo.N_r));  % [Nk, Nr]
        for ll = 1:ch.N_c
            tau_l = (tap_supp(ll)-1)*(1/fq.Bw);
            n_k = (0:fq.N_k-1).';
            ch_u = ch_u + gain(u,ll,np)*exp(-1i*2*pi*tau_l*fq.Bw*n_k/fq.N_k)*Asv(:,u,np).';
        end
                
        u_idx = act_set(u);                       % user id  
        Hf_mtx ( (u_idx-1)*fq.N_k+1:u_idx*fq.N_k,:,np) = ch_u;
        Hf(:,:,u_idx,np) = ch_u;

    end
end

%% channel parameters
ch_para.phi = PHIs;
ch_para.theta = THETAs;
ch_para.tap = tap;
ch_para.gain = gain;
ch_para.Ds = Ds;

end