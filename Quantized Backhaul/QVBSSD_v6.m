function [Rec_sym] = QVBSSD_v6(Y_quant,Y_low,Y_up,H,nvar,M_ord,Data_mod,R_rec,Fo_DFT)
% dealing with quantized obaservations, which only utilize the first-order
% information in VI, can only deal with the time domain transmission
% 2021.4.30 v0
% 2021.5.6 v1 reproduction of the trans
% 2021.5.6 v2 reduce the complexity and rearrange the code
% 2021.5.11 v3 modify the VI to frequency domain
% 2022.7.13 v4 test new,failed
% 2022.7.15 v4 update initialization, success to avoid local minima
% 2022.7.16 v5, modify the prior of the xt to Gaussian distribution - aborted
% 2022.7.20 v6, modify the prior of the xt to time domain 

% reshape dimension, Np is the number of LEOs
[dim1,dim2] = size(Y_quant);
Y_quant = reshape(permute(Y_quant,[2,1]),[dim2,1,dim1]);
Y_low = reshape(permute(Y_low,[2,1]),[dim2,1,dim1]);
Y_up = reshape(permute(Y_up,[2,1]),[dim2,1,dim1]);
R_rec = reshape(permute(R_rec,[2,1]),[dim2,1,dim1]);
Np = length(nvar);  % Np is the number of LEOs
Data_mod = reshape(permute(squeeze(Data_mod(:,1,:)),[2,1]),[size(H,2),1,dim1]);

% parameters predefinition
[Nr,Nu,SCs] = size(H);
Nsym = size(Y_low,2);
% prior parameters
Q = log2(M_ord);            % # bits in each dimension
e_norm = (M_ord==2)+(M_ord~=2)*sqrt((M_ord-1)/6*(2^2));
B_set = unique(qammod((0:M_ord-1)',M_ord,'gray')/e_norm);   % symbol set

% macro definition
Truncdf =@(x,y) normcdf(x)-normcdf(y);
Trunpdf =@(x,y) normpdf(x)-normpdf(y);
% result memory allocation
Rec_sym = zeros(Nu,Nsym,SCs);

% divide the quantization part
central_idx = zeros(Np,1);
central_sel = 1;%randi(Np);
central_idx(central_sel,:) = 1;
Nr_LEO  = Nr/Np;
idx_unquant = ((central_sel-1)*Nr_LEO+1:central_sel*Nr_LEO).'; 
idx_quant = setdiff((1:Nr).',idx_unquant); 


% while mse > 1e-5
for sym = 1:Nsym
    y_low = squeeze(Y_low(:,sym,:));
    y_up = squeeze(Y_up(:,sym,:));
    data_mod = squeeze(Data_mod(:,sym,:));

    % zt = squeeze(R_rec(:,sym,:));
    zt = zeros(Nr,SCs);
    zt(idx_unquant,:)= squeeze(R_rec(idx_unquant,sym,:));
    zt(idx_quant,:)= squeeze(Y_quant(idx_quant,sym,:));  % 即使用准确的zt也会出错
    
    p_xt = 1/length(B_set)*ones(SCs,length(B_set),Nu);

    Xf = zeros(SCs,Nu);  % LMMSE initialization
    for c = 1:SCs
        hall  = H(:,:,c);
        % hall  = H(idx_unquant,:,c);
        Xf(c,:) = ( (hall'*hall+mean(nvar)*eye(Nu))\hall'*zt(:,c) ).';  % right
        % Xf(c,:) = ( (hall'*hall+mean(nvar)*eye(Nu))\hall'*zt(idx_unquant,c) ).';  % right
    end 
    xt_m = Fo_DFT*Xf;  % check this initialization
    xt_var = ones(SCs,Nu); % initialization of variance
    xf_var = diag(mean(xt_var,1));

    % mu_d = zeros(Nr,SCs);
    % for i = 1:iter
    % i = 0; mse = 1; mse_pre = 1;
    % while (mse>1e-3) 
    for i = 1:50
        xtm_pre = xt_m;
        %% module 0-LMMSE estimation of Z
%         for c = 1:SCs
%             hall  = H(:,:,c);        
%             Xf(c,:) = ( (hall'*hall+mean(nvar)*eye(Nu))\hall'*zt(:,c) ).';  % right
%         end  
        for c = 1:SCs  % better than above
            hall  = H(:,:,c);
            % hall  = H(idx_unquant,:,c);
            xf_prior_mean = Xf(c,:).';
            % xf_prior_var =  xf_var;
            xf_prior_var_inv = diag(1./diag(xf_var));
            Xf(c,:) = xf_prior_mean + (xf_prior_var_inv + hall'*(1/mean(nvar))*hall)\hall'*(1/mean(nvar))*(zt(:,c)-hall*xf_prior_mean);
            % Xf(c,:) = xf_prior_mean + (xf_prior_var_inv + hall'*(1/mean(nvar))*hall)\hall'*(1/mean(nvar))*(zt(idx_unquant,c)-hall*xf_prior_mean);
        end

        %% module 1-user wise interference, in time domain

        % QAM  prior
        for u = 1:Nu
            fkt = f10(Fo_DFT',Xf(:,u),p_xt(:,:,u),B_set,mean(nvar),xt_m(:,u));
            % fkt = f10(hall,zt,p_xt,B_set,mean(nvar),xt_m);
            q_xt =  softmax(fkt.').';  % [nUE,nB_set]
            x_m_sc = q_xt*B_set;       % 
            x_v_sc =q_xt*abs(B_set).^2-abs(x_m_sc).^2;
            p_xt(:,:,u) = q_xt;
            xt_m(:,u) = x_m_sc;  % forgotten

            xt_var(:,u) = x_v_sc;  % variance of time domain symbol
        end
        % output of module 1
        Xf = Fo_DFT'*xt_m;
        xf_var = diag(max(mean(xt_var,1),1e-10));
        mu_d = pagemtimes(H,reshape(Xf.',[Nu,1,SCs]));  
        %% module 2
        
        % posterior interference

        for p = 1:Np
            if central_idx(p) == 1
                zt(idx_unquant,:) = zt(idx_unquant,:);
            else
                %                 nvar_np = nvar(p) ;
                nvar_np = mean(nvar);
                nvar_tmp = sqrt(nvar_np/2);
                idx_z = ((p-1)*Nr_LEO+1:p*Nr_LEO).';
                mu_re = real(mu_d(  idx_z,:)); mu_im  = imag(mu_d(  idx_z,:));
                y_low_re = real(y_low(  idx_z,:)); y_low_im =imag(y_low(  idx_z,:));
                y_up_re = real(y_up(  idx_z,:)); y_up_im = imag(y_up(  idx_z,:));

                z_m_re = mu_re + nvar_tmp*Trunpdf( (y_low_re-mu_re)/nvar_tmp, (y_up_re-mu_re)/nvar_tmp )...
                    ./ (Truncdf( (y_up_re-mu_re)/nvar_tmp, (y_low_re-mu_re)/nvar_tmp )+eps);
                z_m_im = mu_im + nvar_tmp*Trunpdf( (y_low_im-mu_im)/nvar_tmp, (y_up_im-mu_im)/nvar_tmp )...
                    ./ (Truncdf( (y_up_im-mu_im)/nvar_tmp, (y_low_im-mu_im)/nvar_tmp )+eps);

                zt(idx_z,:) = z_m_re + 1i*z_m_im;
            end
        end
        % ztt = squeeze(R_rec(:,:,c));
        % zt = ztt;
        mse = norm(xt_m.'-data_mod,2)/norm(data_mod,2);
        % disp(['mse0 = ', num2str(mse), ', i = ', num2str(i),', mse_var = ', num2str(abs((mse-mse_pre)/mse_pre))]);

        % i = i+1;
        if abs(norm(xt_m(:)-xtm_pre(:),2)/norm(xtm_pre(:),2))<1e-20
            break;
        end
        % mse_pre = mse;
        mseall(i) = mse;
    end
    % disp(['mse = ', num2str(mse)]);
    Rec_sym(:,sym,:) = (Fo_DFT'*xt_m).';



end  % symbol end

Rec_sym = squeeze(Rec_sym);

end  % function end

function fkt = f10(hall,zt,p_xt,B_set,nvar,xt_m)
Kall = size(hall,2);
% fkt = zeros(Kall,length(B_set));
% for k = 1:size(hall,2)
%     hk = hall(:,k);
%     term1 = norm(hk,2)^2*abs(B_set).^2;
%     tmp = (conj(xt_m.')*hall'-conj(xt_m(k,:))*hk')*hk;
%     tmp2 = zt'*hk;
%     term2 = 2*real((tmp-tmp2)*B_set);
%     term3 = log(p_xt(k,:)).';
%     fkt(k,:) = (-1/nvar*(term1+term2)+term3).';
% end

term1_ = abs(B_set).^2*sum(abs(hall).^2,1);  % [QAM,Kall]
tmp0_ = repmat(conj(xt_m.')*hall',[Kall,1]);
tmp1_ = diag(conj(xt_m.'))*hall';    % [Kall,]
tmp2_ = (diag((tmp0_-tmp1_)*hall)).';   % [Kall,1]
tmp3_ = zt'*hall;  % [1,Kall]
term2_ = 2*real(B_set*(tmp2_-tmp3_));  % [QAM,Kall]
term3_ = log(p_xt).';  % [QAM,Kall]
fkt = (-1/nvar*(term1_+term2_)+term3_).';


end
