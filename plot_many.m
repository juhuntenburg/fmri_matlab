load('171012.mat')
load('171012all_experiment_stats.mat')
plot(monit.resp_rate)
G1_1=monit.resp_rate(2480:12080); %2660-3440
G1_1t=monit.temp(2480:12080);

load('171013A.mat')
load('171013Aall_experiment_stats.mat')
plot(monit.resp_rate)
G2_1=monit.resp_rate(1790:11390); %1894-2990 4296-4702
G2_1t=monit.temp(1790:11390);

load('171013B.mat')
load('171013Ball_experiment_stats.mat')
plot(monit.resp_rate)
G3_1=monit.resp_rate(2790:11370); 
G3_1t=monit.temp(2790:11370);

load('171013C.mat')
load('171013Call_experiment_stats.mat')
plot(monit.resp_rate)
G3_2=monit.resp_rate(2890:11470); 
G3_2t=monit.temp(2890:11470);

load('171016.mat')
load('171016all_experiment_stats.mat')
plot(monit.resp_rate)
G1_2=monit.resp_rate(2930:11510); %  4695-6110
G1_2t=monit.temp(2930:11510); 

load('171017.mat')
load('171017all_experiment_stats.mat')
plot(monit.resp_rate)
G2_2=monit.resp_rate(2750:11330); 
G2_2t=monit.temp(2750:11330); 

load('171019.mat')
load('171019all_experiment_stats.mat')
plot(monit.resp_rate)
G2_3=monit.resp_rate(2630:11210); 
G2_3t=monit.temp(2630:11210);

load('171020A.mat')
load('171020Aall_experiment_stats.mat')
plot(monit.resp_rate)
G3_3=monit.resp_rate(2932:11512); 
G3_3t=monit.temp(2932:11512);

load('171020B.mat')
load('171020Ball_experiment_stats.mat')
plot(monit.resp_rate)
G1_3=monit.resp_rate(2696:11276); 
G1_3t=monit.temp(2696:11276);

x=1:9601;
y=1021:9601;
plot(x([1:150]),G1_1([1:150]),'Color',[0.0,0.2,1.0]);hold on;
plot(x([962:9601]),G1_1([962:9601]),'Color',[0.0,0.2,1.0]);hold on;
plot(x([1021:1765+1020]),G1_2([1:1765]),'Color',[0.8,0.0,1.0]);hold on;
plot(x([3180+1022:9601]),G1_2([3180+2:9601-1020]),'Color',[0.8,0.0,1.0]);hold on;
plot(x,G1_1,'Color',[0.0,0.2,1.0]);hold on;
plot(y,G1_2,'Color',[0.8,0.0,1.0]);hold on;
plot(y,G1_3,'Color',[0.2,0.8,1.0]);hold on;

plot(x,G2_1,'Color',[1.0,0.0,0.0]);hold on;
plot(y,G2_2,'Color',[1.0,0.8,0.2]);hold on;
plot(y,G2_3,'Color',[0.6,0.0,0.2]);hold on;
plot(y,G3_1,'Color',[0.0,1.0,0.0]);hold on;
plot(y,G3_2,'Color',[0.0,0.4,0.0]);hold on;
plot(y,G3_3,'Color',[0.2,0.6,0.4]);hold on;