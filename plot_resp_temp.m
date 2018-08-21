load('180704C.mat')
%S = sprintf('%s ', monit.resp_rate{:});
%D = sscanf(S, '%f');
plot(monit.resp_rate)
G1_1=monit.resp_rate(280:2608);
G1_1t=monit.temp_new(280:2608);
G1_1e=monit.event(280:2608);

figure(2)
plot(((1:length(G1_1))-1)/60,G1_1,'Color','k');ylim([0,250]);hold on;plot(((1:length(G1_1))-1)/60,250*logical(G1_1e));hold on;
set(gcf,'Position',[38         430        1122         268]);
xlabel('Time [min]');
xlim([0,(length(G1_1)-1)/60])
% line([51 51],get(gca,'YLim'),'Color','r')
% line([16 16],get(gca,'YLim'),'Color','b')
% line([55 55],get(gca,'YLim'),'Color','b')
ylabel('Respiration rate [Resp/min]');
saveas(gcf,'180507Cresp645-11210','tif');

figure(3)
plot(((1:length(G1_1t))-1)/60,G1_1t,'Color','k');hold on;plot(((1:length(G1_1))-1)/60,39*logical(G1_1e));hold on;
set(gcf,'Position',[38         430        1122         268]);
xlabel('Time [min]');
xlim([0,(length(G1_1t)-1)/60])
ylim([32,39])
ylabel('Temperature [\circC]');
% line([51 51],get(gca,'YLim'),'Color','r')
% line([16 16],get(gca,'YLim'),'Color','b')
% line([55 55],get(gca,'YLim'),'Color','b')
saveas(gcf,'180507Ctemp645-11210','tif');