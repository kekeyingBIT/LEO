% ErrorChecking;
% data = bit33;
BerCoded_th = data.BerCoded_th;
BerUncoded_th = data.BerUncoded_th;
Nmse_th = data.Nmse_th;
Acc_th = data.Acc_th;
Pe_th = data.Pe_th;

%%
overhead = [0.036,0.04:0.01:0.12,0.15]*1700;
SNR = -4:4:20;      % snr choose
CE_Sel = char('Perfect-CIR');         % choose JADCE algorithm
DD_Sel = char('Centralized_Quant_LS','Centralized_Quant_LS_mix','Centralized_Quant_VI_3','Centralized_Quant_VI','Centralized_Unquant_LS');                                % choose DD algorithm
CIR_Sel = char('Perfect-CIR');                                  % choose CIR used in DD
DT = char('diversity');                                                    % transmission mode selection
%% plot Ber
startpoint = 2;
BerUncoded = data.BerUncoded_th;
BerCoded = data.BerCoded_th;
codeNum = 1;  % Coded or not
modeNum = 5;  % Centralized or not
CIRNum = 1;   % CIR number
lineSE = zeros(modeNum,CIRNum);
color = char('k','r','b','m','c','g','y');
lineStyle = char('-','-.','--',':');
mark = char('s','p','o','p','v','^');

figure,
for c = 1:codeNum
    if c == 1
        ber = BerUncoded;
        % berperfect = bit55.BerUncoded_th;
    else
        ber = BerCoded;
        % berperfect = bit55.BerCoded_th;
    end

    for k = 1:modeNum


        for j = 1:CIRNum
            lineSE(k,j) = semilogy(overhead(startpoint:end), squeeze(ber(k,startpoint:end)));
            set(lineSE(k,j),'Color',color(k),'LineStyle',lineStyle(c,:),...
                'LineWidth',1.2, 'Marker',mark(j)); hold on; grid on;
        end

    end
end

set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.35)
set(gca,'MinorGridLineStyle',':','MinorGridColor','k',...
    'MinorGridAlpha',0.25)
set(gca,'xcolor','k');
set(gca,'ycolor','k');


LegdStr = char('3 bit backhaul - Alg. 3(MSCTP)',...
    '3 bit backhaul - Alg. 3(MSCBP)',...
    '3 bit backhaul - Alg. 4(MSCTP)',...
    '3 bit backhaul - Alg. 4(MSCBP)',...
    'perfect backhaul - Alg. 3');
legend(LegdStr,'Location','southwest');

xlabel('non-ISI region length');
ylabel('BER');
% axis tight;
