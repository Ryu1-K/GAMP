% close all
LW = 2;
FS = 16;
FN = 'Times New Roman';
MS = 12;
MFC ='w';
CL = 'b';
MK = 'x';

figure(1)
h = semilogy(SIM.EsN0,SIM.BER);
axis([min(SIM.EsN0) max(SIM.EsN0) 10^(-5) 1]);
hold on
grid on
set(h,'LineWidth',LW,'Marker',MK,'MarkerSize',MS,'MarkerFaceColor',MFC)
set(gca,'LineWidth',LW,'FontSize',FS,'FontName',FN,'Ytick',10.^([-5:0]),'PlotBoxAspectRatio',[1,1,1]);
lb(1) = ylabel('BER');
lb(2) = xlabel('\it{E}\rm_s/\it{N}\rm_0 [dB]');
set(lb,'FontSize',FS,'FontName',FN);

saveas(gcf,'./FIG/ber','epsc');
saveas(gcf,'./FIG/ber','png');
saveas(gcf,'./FIG/ber','fig');