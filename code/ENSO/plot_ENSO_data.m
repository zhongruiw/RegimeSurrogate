dt = 0.002; % time step
T = 24000; % total time length, i.e., 1200*2=2400 months=200 years
N = round(T/dt); % total numerical integration steps


[u_3R,h_W_3R,T_C_3R,T_E_3R,tau,I,s1,s2,s3,noise,n,sgm,p,xx] = reference_model(dt,T);

u_lg = 1.5;
h_lg = 150;
t_lg = 7.5;
tau_lg = 5;

dt = 0.002;
k_dt=0.5/dt;
nn=1;

k_dt=0.5/dt;
nn=1;
T_C_3R=smooth(T_C_3R,nn*k_dt);
T_C_3R=T_C_3R(1:k_dt:end);
T_E_3R=smooth(T_E_3R,nn*k_dt);
T_E_3R=T_E_3R(1:k_dt:end);
h_W_3R=smooth(h_W_3R,nn*k_dt);
h_W_3R=h_W_3R(1:k_dt:end);
u_3R=smooth(u_3R,nn*k_dt);
u_3R=u_3R(1:k_dt:end);
I=smooth(I,nn*k_dt);
I=I(1:k_dt:end);
T_C_3R=T_C_3R-mean(T_C_3R);
T_E_3R=T_E_3R-mean(T_E_3R);
h_W_3R=h_W_3R-mean(h_W_3R);
u_3R=u_3R-mean(u_3R);


figure('color','white')
set(gcf,'unit','centimeters','position',[1 1 25 20]);
subplot(2,2,1)
plot(xx,p,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(a) PDF of I','fontsize',16)
xlabel('I')

subplot(2,2,2)
plot(xx,sgm,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(b) Multiplicative noise value \sigma_I(I)','fontsize',16)
xlabel('I')

subplot(2,1,2)
plot(dt:k_dt*dt:N*dt, I,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(c) Time series of I','fontsize',16)
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
xlabel('t')

T_C_3R = T_C_3R * 7.5;
T_E_3R = T_E_3R * 7.5;
h_W_3R = h_W_3R * 150;
u_3R = u_3R * 1.5;
 
disp('Stats of TC model:')
disp([mean(T_C_3R), var(T_C_3R),skewness(T_C_3R),kurtosis(T_C_3R)]);
disp('Stats of TE model:')
disp([mean(T_E_3R), var(T_E_3R),skewness(T_E_3R),kurtosis(T_E_3R)]);



%% observation
data=importdata('./nino3.txt');
nino3=reshape(data(:,2:end)',[],1);
nino3=nino3((1950-1870)*12+1:end-12);

nino3_seasonal = reshape(nino3,12,[]);
nino3_seasonal = var(nino3_seasonal');

data=importdata('./nino4.txt');
nino4=reshape(data(:,2:end)',[],1);
nino4=nino4((1950-1870)*12+1:end-12);

nino4_seasonal = reshape(nino4,12,[]);
nino4_seasonal = var(nino4_seasonal');

data=importdata('./nino34.txt');
nino34=reshape(data(:,2:end)',[],1);
nino34=nino34((1950-1870)*12+1:end-12);

nino34_seasonal = reshape(nino34,12,[]);
nino34_seasonal = var(nino34_seasonal');

disp('Stats of Nino4 obs:')
disp([mean(nino4), var(nino4),skewness(nino4),kurtosis(nino4)]);
disp('Stats of Nino3 obs:')
disp([mean(nino3), var(nino3),skewness(nino3),kurtosis(nino3)]);


T_E_3R_seasonal = reshape(T_E_3R,12,[]);
T_E_3R_seasonal = var(T_E_3R_seasonal');
T_C_3R_seasonal = reshape(T_C_3R,12,[]);
T_C_3R_seasonal = var(T_C_3R_seasonal');


%% spectrum
% T_E_3R_cal=T_E_3R(end-length(nino3)+1:end);
% T_C_3R_cal=T_C_3R(end-length(nino3)+1:end);
% [Period_TE,Density_TE,Frequency_TE]=spectrum_cal(T_E_3R_cal);
% [Period_TC,Density_TC,Frequency_TC]=spectrum_cal(T_C_3R_cal);
% [Period_nino3_obs,Density_nino3_obs,Frequency_nino3_obs]=spectrum_cal(nino3);
% [Period_nino4_obs,Density_nino4_obs,Frequency_nino4_obs]=spectrum_cal(nino4);
% 
% figure('color','white');
% nf=length(Frequency_TE);
% %--------figures
% tnum=1; % tnum=1 year./(Ñù±¾µÄÊ±¼ä¼ä¸ô)
% plot(log10(Period_nino3_obs(2:length(Period_nino3_obs))./tnum),Density_nino3_obs(2:length(Density_nino3_obs)),'k');
% hold on
% plot(log10(Period_nino4_obs(2:length(Period_nino4_obs))./tnum),Density_nino4_obs(2:length(Density_nino4_obs)),'r');
% plot(log10(Period_TE(2:length(Period_TE))./tnum),Density_TE(2:length(Density_TE)),'k--');
% plot(log10(Period_TC(2:length(Period_TC))./tnum),Density_TC(2:length(Density_TC)),'r--');
% % plot(log10(Period(2:nf)./tnum),SK(2:nf),'k','LineWidth',2);
% set(gca,'XTick',[log10(1) log10(12) log10(24) log10(36) log10(48) log10(60) log10(12*7) log10(12*15)])
% set(gca,'XTickLabel',[0 1 2 3 4 5 7 15]);
% set(gca,'xlim',[log10(8) log10(12*20)],'ylim',[0 0.1]);
% % text(log10(Period(33)),Density(33),...
% %     num2str(Period(33)/12),'fontsize',20);
% 
% legend('Nino3\_Obs','Nino4\_Obs','TE','TC');
% %grid off;
% %set(gca,'GridLineStyle','--');
% XX=xlabel('Period (year) ','Fontname','Times','fontsize',14);
% set(gca,'linewidth',.5);
% set(gca,'Fontname','Times','fontsize',14);
% title('The estimate of power spectra of Nino indices','fontsize',20);

figure
subplot(5,4,1:3)
hold on
plot(dt:k_dt*dt:N*dt, T_C_3R ,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('T_C')
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
ylabel('{}^oC')

subplot(5,4,4)
[fi,xx] = ksdensity(T_C_3R);
plot(xx,fi,'g','linewidth',2)
box on
set(gca,'fontsize',12)
subplot(5,4,5:7)
hold on
plot(dt:k_dt*dt:N*dt, T_E_3R ,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('T_E')
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
ylabel('{}^oC')
subplot(5,4,8)
[fi,xx] = ksdensity(T_E_3R);
plot(xx,fi,'r','linewidth',2)
box on
set(gca,'fontsize',12)
subplot(5,4,9:11)
hold on
plot(dt:k_dt*dt:N*dt, h_W_3R ,'b')
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
box on
set(gca,'fontsize',12)
title('h_W')
ylabel('m')
subplot(5,4,12)
[fi,xx] = ksdensity(h_W_3R);
plot(xx,fi,'b','linewidth',2)
box on
set(gca,'fontsize',12)
subplot(5,4,13:15)
hold on
plot(dt:k_dt*dt:N*dt, u_3R ,'k')
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
box on
set(gca,'fontsize',12)
title('u')
ylabel('m/s')
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
xlabel('Time(years)');
subplot(5,4,16)
[fi,xx] = ksdensity(u_3R);
plot(xx,fi,'k','linewidth',2)
box on
set(gca,'fontsize',12)
subplot(5,4,17:19)
hold on
plot(dt:k_dt*dt:N*dt, I ,'k')
set(gca,'xlim',[1 1200]);
set(gca,'xtick',1:10*6:1200,'xticklabel',1:10:200);
box on
set(gca,'fontsize',12)
title('I')
ylabel('')



disp('Stats of hW model:')
disp([mean(h_W_3R), var(h_W_3R),skewness(h_W_3R),kurtosis(h_W_3R)]);
disp('Stats of u model:')
disp([mean(u_3R), var(u_3R),skewness(u_3R),kurtosis(u_3R)]);

load obs_data

disp('Stats of hW obs:')
disp([mean(h_W_obs), var(h_W_obs),skewness(h_W_obs),kurtosis(h_W_obs)]);
disp('Stats of u obs:')
disp([mean(u_obs), var(u_obs),skewness(u_obs),kurtosis(u_obs)]);



spectrum_5


corr1 = corrcoef(T_E_obs,T_C_obs);
disp(['Corr TE TC obs:   ', num2str(corr1(1,2))]);
corr2 = corrcoef(T_E_3R,T_C_3R);
disp(['Corr TE TC model:   ', num2str(corr2(1,2))]);
corr3 = corrcoef(h_W_obs(1:468),u_obs);
disp(['Corr hW u obs:   ', num2str(corr3(1,2))]);
corr4 = corrcoef(h_W_3R,u_3R);
disp(['Corr hW u model:   ', num2str(corr4(1,2))]);
corr5 = corrcoef(T_E_obs(end-480+1:end),h_W_obs);
disp(['Corr TE hW obs:   ', num2str(corr5(1,2))]);
corr6 = corrcoef(T_E_3R,h_W_3R);
disp(['Corr TE hW model:   ', num2str(corr6(1,2))]);
corr7 =  corrcoef(T_C_obs(end-480+1:end-12),u_obs);
disp(['Corr TC u obs:   ', num2str(corr7(1,2))]);
corr8 = corrcoef(T_C_3R,u_3R);
disp(['Corr TC u model:   ', num2str(corr8(1,2))]);