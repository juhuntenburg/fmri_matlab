clear all
close all

folder = '2018_05';
file = '180507C';

resp_rate = load(['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'r_rate.mat']);
r_rate = resp_rate.r_rate;
figure(1)
for i=1:size(r_rate,2)
    subplot(size(r_rate,2),1,i);
    plot(r_rate{i});
    if i==size(r_rate,2)
        xlabel('Time (s)');
    end
    if i==round(size(r_rate,2)/2)
        ylabel('Respiration rate (Resp/m)');
    end
end
saveas(gcf,['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'r_rate_stim'],'tif');

% figure(2)
% for i=1:size(r_rate,2)
%     plot(r_rate{i}); hold on;
% end
% hold off;
% legend('1','2','3','4','5','6','7');

temp = load(['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'tem.mat']);
tem = temp.tem;
figure(3)
for i=1:size(tem,2)
    subplot(size(tem,2),1,i);
    plot(tem{i});
    if i==size(tem,2)
        xlabel('Time (s)');
    end
    if i==round(size(tem,2)/2)
        ylabel('T (\circC)');
    end
end
saveas(gcf,['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'tem_stim'],'tif');

values = [];
for i=1:size(r_rate,2)
    values=[values; r_rate{i}];
    r_rate_stats.mean_resprate(i) = mean(r_rate{i});
    r_rate_stats.std_resprate(i) = std(r_rate{i});
    r_rate_stats.min_resprate(i) = min(r_rate{i});
    r_rate_stats.max_resprate(i) = max(r_rate{i});
end
r_rate_stats.mean_resprate_all = mean(values);
r_rate_stats.std_resprate_all = std(values);
r_rate_stats.min_resprate_all = min(values);
r_rate_stats.max_resprate_all = max(values);
save(['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'r_rate_stats'],'r_rate_stats');

values = [];
for i=1:size(tem,2)
    values=[values; tem{i}];
    tem_stats.mean_tem(i) = mean(tem{i});
    tem_stats.std_tem(i) = std(tem{i});
    tem_stats.min_tem(i) = min(tem{i});
    tem_stats.max_tem(i) = max(tem{i});
end
tem_stats.mean_tem_all = mean(values);
tem_stats.std_tem_all = std(values);
tem_stats.min_tem_all = min(values);
tem_stats.max_tem_all = max(values);
save(['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'tem_stats'],'tem_stats');

monitoring = load(['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file '.mat']);
monit=monitoring.monit;
figure(4);plot(monit.resp_rate); hold on; plot(logical(monit.event)*max(monit.resp_rate));
xlabel('Time (s)');
ylim([min(monit.resp_rate) max(monit.resp_rate)]);
ylabel('Respiration rate (Resp/m)');
saveas(gcf,['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'r_rate'],'tif');

figure(5);plot(monit.temp_new); hold on; plot(logical(monit.event)*max(monit.temp_new));
xlabel('Time (s)');
ylim([min(monit.temp_new) max(monit.temp_new)]);
ylabel('T (\circC)');
saveas(gcf,['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'tem'],'tif');

%%
all_experiment_stats.init_point = 514;
% all_experiment_stats.mid_point1 = 4695;
% all_experiment_stats.mid_point2 = 6110;
% all_experiment_stats.mid_point3 = 4296;
% all_experiment_stats.mid_point4 = 4702;
all_experiment_stats.final_point = 11210;
all_experiment_stats.r_rate_cut=monit.resp_rate(all_experiment_stats.init_point:all_experiment_stats.final_point);
all_experiment_stats.r_rate_mean=mean(monit.resp_rate(all_experiment_stats.init_point:all_experiment_stats.final_point));
all_experiment_stats.r_rate_std=std(monit.resp_rate(all_experiment_stats.init_point:all_experiment_stats.final_point));
all_experiment_stats.tem_cut=monit.temp_new(all_experiment_stats.init_point:all_experiment_stats.final_point);
all_experiment_stats.tem_mean=mean(monit.temp_new(all_experiment_stats.init_point:all_experiment_stats.final_point));
all_experiment_stats.tem_std=std(monit.temp_new(all_experiment_stats.init_point:all_experiment_stats.final_point));

% all_experiment_stats.r_rate_cut=monit.resp_rate([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.final_point]);
% all_experiment_stats.r_rate_mean=mean(monit.resp_rate([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.final_point]));
% all_experiment_stats.r_rate_std=std(monit.resp_rate([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.final_point]));
% all_experiment_stats.tem_cut=monit.temp([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.final_point]);
% all_experiment_stats.tem_mean=mean(monit.temp([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.final_point]));
% all_experiment_stats.tem_std=std(monit.temp([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.final_point]));

% all_experiment_stats.r_rate_cut=monit.resp_rate([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.mid_point3,all_experiment_stats.mid_point4:all_experiment_stats.final_point]);
% all_experiment_stats.r_rate_mean=mean(monit.resp_rate([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.mid_point3,all_experiment_stats.mid_point4:all_experiment_stats.final_point]));
% all_experiment_stats.r_rate_std=std(monit.resp_rate([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.mid_point3,all_experiment_stats.mid_point4:all_experiment_stats.final_point]));
% all_experiment_stats.tem_cut=monit.temp([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.mid_point3,all_experiment_stats.mid_point4:all_experiment_stats.final_point]);
% all_experiment_stats.tem_mean=mean(monit.temp([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.mid_point3,all_experiment_stats.mid_point4:all_experiment_stats.final_point]));
% all_experiment_stats.tem_std=std(monit.temp([all_experiment_stats.init_point:all_experiment_stats.mid_point1,all_experiment_stats.mid_point2:all_experiment_stats.mid_point3,all_experiment_stats.mid_point4:all_experiment_stats.final_point]));
save(['/Users/francisca/Documents/Data/Physiological_Data/' folder '/' file 'all_experiment_stats'],'all_experiment_stats');