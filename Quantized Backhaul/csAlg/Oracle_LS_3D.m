function [X] = Oracle_LS_3D(Y,A,supp,nvar,ce_choice,H_supp)
% 3D oracle-LS 
% supp: [signal-dim] - index mat
[~,N,K] = size(A); % [observation, signal dim, 3rd-dim]
[~,M,~] = size(Y); % [observation, MMV-dim, 3rd-dim];

X = zeros(N,M,K);  % [signal-dim, MMV-dim, 3rd-dim]

for k = 1:K

    A_oracle = A(:,supp,k);
    switch ce_choice
        case 'Oracle-LS'
            x_supp = (A_oracle'*A_oracle+ nvar*eye(length(supp)) )\A_oracle'*Y(:,:,k);
        case 'Oracle-LMMSE'
            x_supp = (A_oracle'*A_oracle+nvar*(H_supp*H_supp'+eps)^(-1))\A_oracle'*Y(:,:,k);
    end
    X(supp,:,k) = x_supp;
end

