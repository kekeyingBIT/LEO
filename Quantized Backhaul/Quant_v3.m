function [Y_mat,Low_mat,Up_mat,agc] = Quant_v3(Y,B)
% 包括自动增益控制AGC
% AGC
tmp_sig = [real(Y(:)),imag(Y(:))];
agc = 1/max(tmp_sig(:));
sig = Y*agc;

% 量化
dim1 = size(sig,1);
dim2=size(sig,2);
dim3=size(sig,3);
QR=zeros(dim1,3,dim2,dim3);
QI=zeros(dim1,3,dim2,dim3);
sig_re = real(sig);
sig_im = imag(sig);
[partition,codebook] = lloyds([sig_re(:);sig_im(:)],2^B);
Upper = +5;
Lower = -5;
for n = 1:dim3
    for  m = 1:dim2
        y=sig(:,m,n);
        [index_re,quants_re,distor_re] = quantiz(real(y),partition,codebook);
        [index_im,quants_im,distor_im] = quantiz(imag(y),partition,codebook);
        threshold  = [Lower partition Upper];
        QR(:,1,m,n) = quants_re.';
        QR(:,2,m,n) = threshold(index_re+1);
        QR(:,3,m,n) = threshold(index_re+2);
        QI(:,1,m,n) = quants_im.';
        QI(:,2,m,n) = threshold(index_im+1);
        QI(:,3,m,n) = threshold(index_im+2);
    end
end
 bb_sig = squeeze(QR(:,1,:,:)+1i*QI(:,1,:,:));
 bb_sig_low = squeeze(QR(:,2,:,:)+1i*QI(:,2,:,:));
 bb_sig_up = squeeze(QR(:,3,:,:)+1i*QI(:,3,:,:));
 Y_mat = bb_sig/agc;  % time domain signal
 Low_mat = bb_sig_low/agc;
 Up_mat = bb_sig_up/agc;
