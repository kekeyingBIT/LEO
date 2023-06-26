% ErrorChecking;

BerCoded_th = data.BerCoded_th;
BerUncoded_th = data.BerUncoded_th;
Nmse_th = data.Nmse_th;
Acc_th = data.Acc_th;
Acc_th2 = data.Acc_th2;
Pe_th = data.Pe_th;
Pe_th2 = data.Pe_th2;
%%
overhead = [0.15:0.05:0.45]*100; % leo.NpAll = 1:7;
CE_Sel = char('SOMP','Oracle-LS','OAMP-1','OAMP-2');         % choose JADCE algorithm
DD_Sel = char('Distributed','Centralized');                                % choose DD algorithm
CIR_Sel = char('SOMP','Oracle-LS','OAMP-1','OAMP-2');                                                  % choose CIR used in DD
DT = char('diversity');                                                    % transmission mode selection
Coded_Sel = char('No');                                                  % Channel coding scheme, No/LDPC

%% define figure parameter
color = char('k','r','b','m','c','g','y');
lineStyle = char(':','-.','--','-','-');
mark = char('s','p','o','p','v','^');

%% plot nmse
modeNum = 4;
modeIdx = [1,2,3,4];
% gIdx = [1,2,4:1:11];
gIdx = [1:7];
startpoint = 1;
lineSE = zeros(modeNum,1);
for k = 1:modeNum
    figure(1),
    nmse = data.Nmse_th;
    lineSE(k) = plot(overhead(gIdx), nmse(modeIdx(k),gIdx));
    set(lineSE(k),'Color',color(k),'LineStyle',lineStyle(k,:),...
        'LineWidth',1.2,'Marker',mark(k));
    hold on; grid on;
end
figure(1)
set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.35)
set(gca,'MinorGridLineStyle',':','MinorGridColor','k',...
    'MinorGridAlpha',0.25)
set(gca,'xcolor','k');
set(gca,'ycolor','k');
% LegdStr = char('SOMP','Oracle-LS','Oracle-LMMSE','OAMP','OAMP-2');
LegdStr = char('SOMP','Oracle-LS','Alg. 1 ','Alg. 1 + Alg. 2');
legend(LegdStr,'Location','southeast');
xlabel('Number of active UTs');
ylabel('NMSE (dB)');
% axis tight;
%% plot aud/Pe

%% plot Ber
startpoint = 1;
BerUncoded = data.BerUncoded_th;
BerCoded = data.BerCoded_th;
codeNum = 1;
modeNum = 2;
CIRNum = 4;
lineSE = zeros(modeNum,CIRNum);
color = char('k','r','b','m','c','g','y');
lineStyle = char('-.','--','-',':');
mark = char('s','p','o','p','v','^');

for c = 1:codeNum
    if c == 1
        ber = BerUncoded;
    else
        ber = BerCoded;
    end

    for k = 1:modeNum
        figure(2),
        
        for j = 1:CIRNum
            lineSE(k,j) = semilogy(overhead(gIdx), squeeze(ber(k,j,gIdx)));
            set(lineSE(k,j),'Color',color(j),'LineStyle',lineStyle(k,:),...
                'LineWidth',1.2, 'Marker',mark(j)); hold on; grid on;
        end

    end
end
figure(2)
set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.35)
set(gca,'MinorGridLineStyle',':','MinorGridColor','k',...
    'MinorGridAlpha',0.25)
set(gca,'xcolor','k');
set(gca,'ycolor','k');


LegdStr = char('SOMP',...
    'Oracle-LS',...
    'Alg. 1',...
    'Alg. 1 + Alg. 2');
legend(LegdStr,'Location','southeast');

xlabel('Number of active UTs');
ylabel('BER');
axis tight;

%% plot Pe
figure(3)
csi_idx = [1,3,4];
lineSE = zeros(length(csi_idx),2);
for j = 1:length(csi_idx)
    lineSE(j,2) = semilogy(overhead(startpoint:end), squeeze(Pe_th(csi_idx(j),:)));
    set(lineSE(j,2),'Color',color(j),'LineStyle',lineStyle(2,:),...
        'LineWidth',1.2, 'Marker',mark(j)); hold on; grid on;
end

for j = 1:length(csi_idx)
    lineSE(j,1) = semilogy(overhead(startpoint:end), squeeze(Pe_th2(csi_idx(j),:)));
    set(lineSE(j,1),'Color',color(j),'LineStyle',lineStyle(1,:),...
        'LineWidth',1.2, 'Marker',mark(j)); hold on; grid on;
end

set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.35)
set(gca,'MinorGridLineStyle',':','MinorGridColor','k',...
    'MinorGridAlpha',0.25)
set(gca,'xcolor','k');
set(gca,'ycolor','k');

LegdStr = char('SOMP',...
    'Alg. 1',...
    'Alg. 1 + Alg. 2');
legend(LegdStr,'Location','southwest');


xlabel('Number of active UTs');
ylabel('AEP');
axis tight;
