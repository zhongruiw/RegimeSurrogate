%% PDFs_Spectrums
% figure('color','white')
% set(gcf,'unit','centimeters','position',[10 5 18 10])
% hold on
% Fs = 12;
% L = length(nino3);%L_temp = L;
% sm=3;
% NFFT = 2^nextpow2(L);NFFT_temp = NFFT;
% Y_u1 = fft(nino3,NFFT)/L; Y_u1 = smooth(Y_u1,sm);  Y_u1 = Y_u1/norm(Y_u1(10:end/2));
% f_u1 = Fs/2*linspace(0,1,NFFT/2+1);
% L = length(u_3R);
% NFFT = NFFT_temp;2^nextpow2(L);f_u2 = Fs/2*linspace(0,1,NFFT/2+1);
% gap = (f_u1(2)-f_u1(1))/(f_u2(2)-f_u2(1));
% Y_u2 = fft(u_3R,NFFT)/L; Y_u2 = smooth(Y_u2,sm);  Y_u2 = Y_u2/norm(Y_u2(10*gap:gap:NFFT/2));
% plot(f_u2,2*abs(Y_u2(1:NFFT/2+1)),'b','linewidth',1);
% plot(f_u1,2*abs(Y_u1(1:NFFT_temp/2+1)),'r','linewidth',2);
% set(gca,'xscale','log')
% box on
% set(gca,'fontsize',9)
% xlim([0.2,2])
% h=legend('Model','Obs','location','northwest','Orientation','vertical');
% set(h,'box','off');
% title('Spectrum of Tau','fontsize',10)
% xlabel('Year')
% set(gca,'xtick',[0.1,0.2,1/3,0.5,1,2])
% set(gca,'xticklabel',[10,5,3,2,1,0.5])
% set(gca,'ytick',0:0.1:0.5);

figure
subplot(2,4,1)
ts = 1*12:1*12:480*12;
Fs = 1*12;
L = length(ts);
NFFT = 2^nextpow2(L); 
Y_y1 = fft(T_E_3R(1:480),NFFT)/L;
f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
tpp_y2 = 2*abs(Y_y1(1:NFFT/2+1));
hold on
plot(1./f_y1(end:-1:1) ,tpp_y2(end:-1:1),'b','linewidth',2) 
set(gca,'xscale','log')
xlim([0.8,10]);
set(gca,'fontsize',12)
box on
xlabel('Year')
set(gca,'xTick',[1:6,8,10])
set(gca,'xTicklabel',  [1:6,8,10]);
hold on
ts = 1*12:1*12:480*12;
Fs = 1*12;
L = length(ts);
NFFT = 2^nextpow2(L); 
Y_y1 = fft(nino3,NFFT)/L;
f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
plot(1./f_y1(end:-1:1),tpp_y(end:-1:1),'r','linewidth',2);
temp_y = 2*abs(Y_y1(1:NFFT/2+1));
tp_y = temp_y; 
title('Nino 3')



figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 10])
subplot(2,4,1)
hold on
Fs = 12;
L = length(nino3);%L_temp = L;
sm=3;
NFFT = 2^nextpow2(L);NFFT_temp = NFFT;
Y_u1 = fft(nino3,NFFT)/L; Y_u1 = smooth(Y_u1,sm);  Y_u1 = Y_u1/norm(Y_u1(10:end/2));
f_u1 = Fs/2*linspace(0,1,NFFT/2+1);
L = length(T_E_3R);
NFFT = NFFT_temp;2^nextpow2(L);f_u2 = Fs/2*linspace(0,1,NFFT/2+1);
gap = (f_u1(2)-f_u1(1))/(f_u2(2)-f_u2(1));
Y_u2 = fft(T_E_3R,NFFT)/L; Y_u2 = smooth(Y_u2,sm);  Y_u2 = Y_u2/norm(Y_u2(10*gap:gap:NFFT/2));
plot(f_u2,2*abs(Y_u2(1:NFFT/2+1)),'b','linewidth',1);
plot(f_u1,2*abs(Y_u1(1:NFFT_temp/2+1)),'r','linewidth',2);
set(gca,'xscale','log')
box on
set(gca,'fontsize',9)
xlim([0.2,2])
h=legend('Model','Obs','location','northwest','Orientation','vertical');
set(h,'box','off');
title('(a) Spectrum of T_E','fontsize',10)
xlabel('Year')
set(gca,'xtick',[0.1,0.2,1/3,0.5,1,2])
set(gca,'xticklabel',[10,5,3,2,1,0.5])
set(gca,'ytick',0:0.1:0.5);

subplot(2,4,5)
hold on
Fs = 12;
L = length(nino4);%L_temp = L;
NFFT = 2^nextpow2(L);NFFT_temp = NFFT;
Y_u1 = fft(nino4,NFFT)/L; Y_u1 = smooth(Y_u1,sm);  Y_u1 = Y_u1/norm(Y_u1(10:end/2));
f_u1 = Fs/2*linspace(0,1,NFFT/2+1);
L = length(T_C_3R);
NFFT = NFFT_temp;2^nextpow2(L);
gap = (f_u1(2)-f_u1(1))/(f_u2(2)-f_u2(1));
Y_u2 = fft(T_C_3R,NFFT)/L; Y_u2 = smooth(Y_u2,sm);  Y_u2 = Y_u2/norm(Y_u2(10*gap:gap:NFFT/2));
plot(f_u2,2*abs(Y_u2(1:NFFT/2+1)),'b','linewidth',1);
plot(f_u1,2*abs(Y_u1(1:NFFT_temp/2+1)),'r','linewidth',2);
set(gca,'xscale','log')
box on
set(gca,'fontsize',9)
xlim([0.2,2])
title('(b) Spectrum of T_C','fontsize',10)
% xlabel('Year')
set(gca,'xtick',[0.1,0.2,1/3,0.5,1,2])
set(gca,'xticklabel',[10,5,3,2,1,0.5])
set(gca,'ytick',0:0.1:0.5);


ACF_T_E_model = autocorr(T_E_3R,60);
ACF_T_C_model = autocorr(T_C_3R,60);
ACF_T_E_obs = autocorr(nino3,60);
ACF_T_C_obs = autocorr(nino4,60);



[PDF_T_E_obs, xx_T_E] = ksdensity(nino3);
[PDF_T_E_model, xx_T_E] = ksdensity(T_E_3R,xx_T_E);
[PDF_T_C_obs, xx_T_C] = ksdensity(nino4);
[PDF_T_C_model, xx_T_C] = ksdensity(T_C_3R,xx_T_C);
subplot(2,4,2)
hold on
plot(xx_T_E,PDF_T_E_obs,'r','linewidth',2)
plot(xx_T_E,PDF_T_E_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(c) PDF of T_E','fontsize',10)
% xlabel('¡ãC')
set(gca,'xtick',-2:2:4);
set(gca,'ytick',0:0.1:0.5);

subplot(2,4,6)
hold on
plot(xx_T_C,PDF_T_C_obs,'r','linewidth',2)
plot(xx_T_C,PDF_T_C_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(d) PDF of T_C','fontsize',10)
%xlabel('¡ãC')
set(gca,'xtick',-2:2:4);
set(gca,'ytick',0:0.2:0.6);

%%%%%%






load obs_data
[PDF_h_W_obs, xx_h_W] = ksdensity(h_W_obs);
[PDF_h_W_model, xx_h_W] = ksdensity(h_W_3R,xx_h_W);
[PDF_u_obs, xx_u] = ksdensity(u_obs);
[PDF_u_model, xx_u] = ksdensity(u_3R,xx_u);
ACF_h_W_obs = autocorr(h_W_obs,60);
ACF_h_W_model = autocorr(h_W_3R,60);
ACF_u_obs = autocorr(u_obs,60);
ACF_u_model = autocorr(u_3R,60);
% figure
subplot(2,4,4)
hold on
plot(xx_h_W,PDF_h_W_obs,'r','linewidth',2)
plot(xx_h_W,PDF_h_W_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(g) PDF of h_W','fontsize',10)
xlabel('m')
set(gca,'xtick',-40:40:40);
set(gca,'ytick',0:0.01:0.03);

subplot(2,4,8)
hold on
plot(xx_u,PDF_u_obs,'r','linewidth',2)
plot(xx_u,PDF_u_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(h) PDF of u','fontsize',10)
xlabel('m/s')
set(gca,'xtick',-0.4:0.4:0.4);
set(gca,'ytick',0:3);


subplot(2,4,3)
hold on
plot(1:12,nino3_seasonal,'r','linewidth',2)
plot(1:12,T_E_3R_seasonal,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(e) Variance of T_E','fontsize',10)
ylabel('¡ãC')
% xlabel('Calendar month')
set(gca,'xtick',1:2:12,'xticklabel',{'J','M','M','J','S','N'})
set(gca,'xlim',[1,12],'ylim',[0.3 1.5]);
set(gca,'ytick',-0.8:0.5:1.6);
% set(gca,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});

subplot(2,4,7)
hold on
plot(1:12,nino4_seasonal,'r','linewidth',2)
plot(1:12,T_C_3R_seasonal,'b','linewidth',2)
box on
set(gca,'fontsize',9)
ylabel('¡ãC')
title('(f) Variance of T_C','fontsize',10)
xlabel('Calendar month')
set(gca,'xlim',[1,12],'ylim',[0.2 0.7]);
set(gca,'xtick',1:2:12,'xticklabel',{'J','M','M','J','S','N'})
set(gca,'ytick',0.2:0.1:0.6);

%% spectrum for Nino3.4 index

load reg_te_tc_hw_ssta
load lon_sst

Hov=zeros(length(lon_sst),length(T_E_3R));
xx=[ones(size(T_E_3R)),T_E_3R,T_C_3R];
for i=1:length(lon_sst)
    Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
end
n34_model=squeeze(nanmean(Hov(72:121,:),1));

figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 8])
subplot(1,3,1)
hold on
Fs = 12;
L = length(nino34);%L_temp = L;
sm=3;
NFFT = 2^nextpow2(L);NFFT_temp = NFFT;
Y_u1 = fft(nino34,NFFT)/L; Y_u1 = smooth(Y_u1,sm);  Y_u1 = Y_u1/norm(Y_u1(10:end/2));
f_u1 = Fs/2*linspace(0,1,NFFT/2+1);
L = length(n34_model);
NFFT = NFFT_temp;2^nextpow2(L);f_u2 = Fs/2*linspace(0,1,NFFT/2+1);
gap = (f_u1(2)-f_u1(1))/(f_u2(2)-f_u2(1));
Y_u2 = fft(n34_model,NFFT)/L; Y_u2 = smooth(Y_u2,sm);  Y_u2 = Y_u2/norm(Y_u2(10*gap:gap:NFFT/2));
plot(f_u2,2*abs(Y_u2(1:NFFT/2+1)),'b','linewidth',1);
plot(f_u1,2*abs(Y_u1(1:NFFT_temp/2+1)),'r','linewidth',2);
set(gca,'xscale','log')
box on
set(gca,'fontsize',9)
xlim([0.2,2])
h=legend('Model','Obs','location','northwest','Orientation','vertical');
set(h,'box','off');
title('(a) Spectrum of Nino3.4 index','fontsize',10)
xlabel('Year')
set(gca,'xtick',[0.1,0.2,1/3,0.5,1,2])
set(gca,'xticklabel',[10,5,3,2,1,0.5])
set(gca,'ytick',0:0.1:0.5);

[PDF_n34_obs, xx_n34] = ksdensity(nino34);
[PDF_n34_model, xx_n34] = ksdensity(n34_model,xx_n34);
subplot(1,3,2)
hold on
plot(xx_n34,PDF_n34_obs,'r','linewidth',2)
plot(xx_n34,PDF_n34_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(b) PDF of Nino3.4 index','fontsize',10)
% xlabel('¡ãC')
set(gca,'xtick',-2:2:4);
set(gca,'ytick',0:0.1:0.5);

n34_model_seasonal = reshape(n34_model,12,[]);
n34_model_seasonal = var(n34_model_seasonal');

subplot(1,3,3)
hold on
plot(1:12,nino34_seasonal,'r','linewidth',2)
plot(1:12,n34_model_seasonal,'b','linewidth',2)
box on
set(gca,'fontsize',9)
title('(c) Variance of Nino3.4 index','fontsize',10)
ylabel('¡ãC')
% xlabel('Calendar month')
set(gca,'xtick',1:2:12,'xticklabel',{'J','M','M','J','S','N'})
set(gca,'xlim',[1,12],'ylim',[0.3 1.5]);
set(gca,'ytick',-0.8:0.5:1.6);


figure

subplot(2,2,1)
hold on
plot([0:60]/12, ACF_T_E_obs,'r','linewidth',2)
plot([0:60]/12, ACF_T_E_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
legend('Nino 3','T_E')
title('ACF','fontsize',9)
subplot(2,2,3)
hold on
plot([0:60]/12, ACF_T_C_obs,'r','linewidth',2)
plot([0:60]/12, ACF_T_C_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
legend('Nino 4','T_C')
subplot(2,2,2)
hold on
plot(0:60,ACF_h_W_obs,'r','linewidth',2)
plot(0:60,ACF_h_W_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
legend('hW Obs','hW model')
title('ACF h_W')
subplot(2,2,4)
hold on
plot(0:60,ACF_u_obs,'r','linewidth',2)
plot(0:60,ACF_u_model,'b','linewidth',2)
box on
set(gca,'fontsize',9)
legend('u Obs','u model')
title('ACF u')

%% Model_Obs_TimeSeries
h_W_obs_new = h_W_obs(1:end-12);
u_obs_new = u_obs;
T_E_obs_new = T_E_obs(30*12+1:end-12);
T_C_obs_new = T_C_obs(30*12+1:end-12);
t_obs = 1980+1/12:1/12:2019;
LL = length(t_obs);
range_model = 1+12*220:LL+12*220;
% range_model = 1+12*480:LL+12*480;
% range_model = 1+12*570:LL+12*570;
range_tau = range_model(1)*k_dt:range_model(end)*k_dt;
t_model = range_model/12;
tau_model = range_tau/k_dt/12;

figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 15])
subplot(5,1,1)
hold on
yyaxis left
plot(t_obs,T_E_obs_new,'r','linewidth',2)
set(gca,'ycolor','r')
xlim([t_obs(1),t_obs(end)])
ylim([-3.5,3.5])
ylabel('^oC')
yyaxis right
plot(t_obs,T_C_obs_new,'g','linewidth',2)
set(gca,'ycolor','g')
ylim([-3.5,3.5])
h=legend('T_E','T_C','orientation','horizontal','location','southwest');
set(h,'box','off');
xlim([t_obs(1),t_obs(end)])
set(gca,'fontsize',9)
ylabel('^oC')
box on
text(1973.5,0,'(a) Obs','fontsize',9)
grid on
grid(gca,'minor')
title('Comparison of the observational time series and model simulations','fontsize',12)

subplot(5,1,2)
hold on
yyaxis left
plot(t_obs,h_W_obs_new,'b','linewidth',2)
set(gca,'ycolor','b')
xlim([t_obs(1),t_obs(end)])
ylim([-50,50])
ylabel('m')
yyaxis right
plot(t_obs,u_obs_new,'k','linewidth',2)
set(gca,'ycolor','k')
ylim([-0.5,0.5])
ylabel('m/s')
h=legend('h_W','u','orientation','horizontal','location','northwest');
set(h,'box','off');
xlim([t_obs(1),t_obs(end)])
set(gca,'fontsize',9)
box on
text(1973.5,0,'(b) Obs','fontsize',9)
grid on
grid(gca,'minor')

subplot(5,1,3)
hold on
yyaxis left
plot(t_model,T_E_3R(range_model),'r','linewidth',2)
set(gca,'ycolor','r')
xlim([t_model(1),t_model(end)])
ylim([-3.5,3.5])
ylabel('^oC')
yyaxis right
plot(t_model,T_C_3R(range_model),'g','linewidth',2)
set(gca,'ycolor','g')
ylim([-3.5,3.5])
ylabel('^oC')
h=legend('T_E','T_C','orientation','horizontal','location','southwest');
set(h,'box','off');
xlim([t_model(1),t_model(end)])
set(gca,'fontsize',9)
box on
text(t_model(1)-6.5,0,'(c) Model','fontsize',9)
grid on
grid(gca,'minor')

subplot(5,1,4)
hold on
yyaxis left
plot(t_model,h_W_3R(range_model),'b','linewidth',2)
set(gca,'ycolor','b')
xlim([t_model(1),t_model(end)])
ylim([-50,50])
ylabel('m')
yyaxis right
plot(t_model,u_3R(range_model),'k','linewidth',2)
set(gca,'ycolor','k')
ylim([-0.5,0.5])
ylabel('m/s')
h=legend('h_W','u','orientation','horizontal','location','northwest');
set(h,'box','off');
xlim([t_model(1),t_model(end)])
set(gca,'fontsize',9)
box on
text(t_model(1)-6.5,0,'(d) Model','fontsize',9)
grid on
grid(gca,'minor')

subplot(5,1,5)
hold on
yyaxis left
plot(tau_model,tau(range_tau)*5,'c','linewidth',2)
xlim([t_model(1),t_model(end)])
ylim([-15,15])
set(gca,'ycolor','c')
ylabel('m/s')
yyaxis right
plot(t_model,I(range_model),'m','linewidth',2)
xlim([t_model(1),t_model(end)])
set(gca,'ycolor','m')
ylim([0,4])
% ylabel('^oC')
h=legend('\tau','I','orientation','horizontal','location','southwest');
set(h,'box','off');
box on
set(gca,'fontsize',9)
text(t_model(1)-6.5,2,'(e) Model','fontsize',9)
grid on
grid(gca,'minor')
xlabel('t')


%% Model_hov_particular_events
load reg_te_tc_hw_ssta
load lon_sst



year_add = 20;
t_temp = 1980+1/12:1/12:1980+year_add;%20 years
LL = length(t_temp);

window = 12;
total_loop = (year_add-1) * 12/window;
figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 15])
colormap(jet)
for j = 1:5
    if j == 1
        range_model = 1+12*220:LL+12*220;
%                 range_model = 1+12*900:LL+12*900;
    elseif j == 2
        range_model = 1+12*240:LL+12*240;
%                 range_model = 1+12*920:LL+12*920;
    elseif j == 3
        range_model = 1+12*1410:LL+12*1410;
%                 range_model = 1+12*940:LL+12*940;
    elseif j == 4
        range_model = 1+12*750:LL+12*750;
%                 range_model = 1+12*960:LL+12*960;
    elseif j == 5
        range_model = 1+12*900:LL+12*900;
    end
    range_tau = range_model(1)*k_dt:range_model(end)*k_dt;
    t_model = range_model/12;
    tau_model = range_tau/k_dt/12;
    Hov=zeros(length(lon_sst),length(range_model));
    
    xx=[ones(size(T_E_3R(range_model))),T_E_3R(range_model),T_C_3R(range_model)];
    for i=1:length(lon_sst)
        Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
    end
    subplot(1,5,j)
    [xx,yy] = meshgrid(t_model,lon_sst);
    contourf(yy,xx,Hov,30,'linestyle','none')
    hold on
    plot([180 180],[t_model(1) t_model(end)],'m--','linewidth',2);
    temp_tau = range_tau(1:k_dt/3:end);
    plot(180+tau(temp_tau)*20,tau_model(1:k_dt/3:end),'k','linewidth',0.5);
    for k = 1:total_loop
        if mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 1.0 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
            plot([120,120],[t_model(1-4+k*window),t_model(1+1+window*k)],'r','linewidth',10)
        elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
            plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'m','linewidth',10)
        elseif mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
            plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'color',[255 97 0]/255,'linewidth',10)
        elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5 || mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5
            plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'b','linewidth',10)
        end
    end
    %     colorbar
    set(gca,'xlim',[120 280]);
    set(gca,'xtick',120:60:280);
    if j == 3
        set(gca,'yticklabel',402:2:420);
    end
    if j == 4
        set(gca,'yticklabel',442:2:460);
    end
    if j == 5
        set(gca,'yticklabel',502:2:520);
    end
    
    xlabel('longitude');
    set(gca,'fontsize',9,'linewidth',1);
    %     title('Evolutions of SSTA and wind burst','fontsize',9);
    caxis([-3,3])
    %     colorbar('eastoutside','position',[.92,0.12,0.005,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');
    if j == 1
        ylabel('model year');
        %         text(460,241,'Hovmoller diagram of the standard run ','fontsize',16)
    end
end
sgtitle('Hovmoller diagrams of the standard run','fontsize',12);
colorbar('eastoutside','position',[.92,0.12,0.01,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');


% %% multi-year El Nino and La Nina
range_model = 1:length(T_E_3R);
Hov=zeros(length(lon_sst),length(range_model));

xx=[ones(size(T_E_3R(range_model))),T_E_3R(range_model),T_C_3R(range_model)];
for i=1:length(lon_sst)
    Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
end


total_loop=length(T_E_3R)/12-1;
sst_max_long=zeros(total_loop,2);
sst_min_long_LN=zeros(total_loop,2);
EN=zeros(total_loop,1);
EPEN=zeros(total_loop,1);
CPEN=zeros(total_loop,1);
LN=zeros(total_loop,1);
EEN=zeros(total_loop,1);%extreme El Nino. N3>2.
for k = 1:total_loop
    if max(T_E_3R(-8+k*window:+3+k*window))>=2.5
        EEN(k)=1;
    end
    if mean(T_E_3R(-0+k*window:+2+k*window)) > 0.5 || mean(T_C_3R(-0+k*window:+2+k*window)) > 0.5
        EN(k)=1;
        Hov_cal=mean(Hov(:,-0+k*window:+2+k*window),2);
        sst_max_long(k,1)=max(Hov_cal);
        sst_max_long(k,2)=find(Hov_cal==max(Hov_cal));
        if mean(T_E_3R(-0+k*window:+2+k*window))>mean(T_C_3R(-0+k*window:+2+k*window))
            EPEN(k)=1;
        else
            CPEN(k)=1;
        end
    elseif mean(T_E_3R(-0+k*window:+2+k*window)) < -0.5 || mean(T_C_3R(-0+k*window:+2+k*window)) < -0.5
        LN(k)=1;  
        Hov_cal=mean(Hov(:,-0+k*window:+2+k*window),2);
        sst_min_long_LN(k,1)=min(Hov_cal);
        sst_min_long_LN(k,2)=find(Hov_cal==min(Hov_cal));
    end
end
MYEN=EN(1:end-1).*EN(2:end);
MYLN=LN(1:end-1).*LN(2:end);
for i=1:length(MYEN)-1
    if MYEN(i)==1 && MYEN(i+1)==1
        MYEN(i+1)=0;
    end
    if MYLN(i)==1 && MYLN(i+1)==1
        MYLN(i+1)=0;
    end
end

sum(EN)
sum(LN)
sum(EPEN)
sum(CPEN)
sum(MYEN)
sum(MYLN)
sum(EEN)

sst_max_long_num=zeros(length(lon_sst),1);
for i=1:length(lon_sst)
    sst_max_long_num(i)=sum(sst_max_long(:,2)==i);
end
sst_max_num=zeros(50,1);
for i=1:50
    for k=1:total_loop
        if sst_max_long(k,1)>(i-1)*0.1 && sst_max_long(k,1)<=i*0.1
            sst_max_num(i)=sst_max_num(i)+1;
        end
    end
end

sst_max_long_num_LN=zeros(length(lon_sst),1);
for i=1:length(lon_sst)
    sst_max_long_num_LN(i)=sum(sst_min_long_LN(:,2)==i);
end
sst_min_num_LN=zeros(50,1);
for i=1:50
    for k=1:total_loop
        if sst_min_long_LN(k,1)<(i-1)*(-0.1) && sst_min_long_LN(k,1)>=i*(-0.1)
            sst_min_num_LN(i)=sst_min_num_LN(i)+1;
        end
    end
end

figure('color','white')
set(gcf,'unit','centimeters','position',[20 5 18 15])
sgtitle({'Bivariate distribution of DJF El Nino SSTA peaks';'2000yr Ctrl, averaged 5S-5N'},'fontsize',14);
subplot(3,3,1:2)
plot(lon_sst,sst_max_long_num,'k','linewidth',2);
hold on
plot([160 160],[0 100],'k--','linewidth',0.5);
plot([210 210],[0 100],'k--','linewidth',0.5);
plot([270 270],[0 100],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[0 100]);
set(gca,'xtick',120:20:280,'xticklabel',[]);
set(gca,'ytick',0:20:100);
ylabel('event count');
set(gca,'fontsize',14);


subplot(3,3,[4:5 7:8])
scatter(lon_sst(sst_max_long(find(sst_max_long(:,2)~=0),2)),sst_max_long(find(sst_max_long(:,2)~=0),1),'k*');
hold on
plot([160 160],[0 6],'k--','linewidth',0.5);
plot([210 210],[0 6],'k--','linewidth',0.5);
plot([270 270],[0 6],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[0 5]);
set(gca,'xtick',120:20:280);
set(gca,'ytick',0:1:5);
ylabel('peak SSTA (degC)');
xlabel('longitude of peak SSTA');
box on
set(gca,'fontsize',14);


subplot(3,3,[6 9])
% xx=0:0.1:4.9;
plot(sst_max_num,0:0.1:4.9,'k','linewidth',2);
set(gca,'xlim',[0 100],'ylim',[0 5]);
set(gca,'xtick',0:20:100);
set(gca,'ytick',0:1:5);
xlabel('event count');
set(gca,'fontsize',14);
title([num2str(sum(EN)) ' total warm events'],'fontsize',10);


figure('color','white')
set(gcf,'unit','centimeters','position',[20 5 18 15])
sgtitle({'Bivariate distribution of DJF La Nina SSTA peaks';'2000yr Ctrl, averaged 5S-5N'},'fontsize',14);
subplot(3,3,1:2)
plot(lon_sst,sst_max_long_num_LN,'k','linewidth',2);
hold on
plot([160 160],[0 200],'k--','linewidth',0.5);
plot([210 210],[0 200],'k--','linewidth',0.5);
plot([270 270],[0 200],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[0 200]);
set(gca,'xtick',120:20:280,'xticklabel',[]);
set(gca,'ytick',0:50:200);
ylabel('event count');
set(gca,'fontsize',14);


subplot(3,3,[4:5 7:8])
scatter(lon_sst(sst_min_long_LN(find(sst_min_long_LN(:,2)~=0),2)),sst_min_long_LN(find(sst_min_long_LN(:,2)~=0),1),'k*');
hold on
plot([160 160],[-5 0],'k--','linewidth',0.5);
plot([210 210],[-5 0],'k--','linewidth',0.5);
plot([270 270],[-5 0],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[-5 0]);
set(gca,'xtick',120:20:280);
set(gca,'ytick',-5:1:0);
ylabel('peak SSTA (degC)');
xlabel('longitude of peak SSTA');
box on
set(gca,'fontsize',14);


subplot(3,3,[6 9])
% xx=0:0.1:4.9;
plot(sst_min_num_LN,0:-0.1:-4.9,'k','linewidth',2);
set(gca,'xlim',[0 100],'ylim',[-5 0]);
set(gca,'xtick',0:20:100);
set(gca,'ytick',-5:1:0);
xlabel('event count');
set(gca,'fontsize',14);
title([num2str(sum(LN)) ' total cold events'],'fontsize',10);

% ylabel('event count');
% T_E_3R_seasonal = reshape(T_E_3R(range_model),12,[]);
% T_E_3R_seasonal = var(T_E_3R_seasonal');
% T_C_3R_seasonal = reshape(T_C_3R(range_model),12,[]);
% T_C_3R_seasonal = var(T_C_3R_seasonal');
%
% figure
% subplot(2,1,1)
% hold on
% plot(1:12,nino3_seasonal,'r','linewidth',2)
% plot(1:12,T_E_3R_seasonal,'r--','linewidth',2)
% box on
% set(gca,'fontsize',9)
% subplot(2,1,2)
% hold on
% plot(1:12,nino4_seasonal,'g','linewidth',2)
% plot(1:12,T_C_3R_seasonal,'g--','linewidth',2)
% box on
% set(gca,'fontsize',9)

% disp('Stats of TC model:')
% disp([mean(T_C_3R(range_model)), var(T_C_3R(range_model)),skewness(T_C_3R(range_model)),kurtosis(T_C_3R(range_model))]);
% disp('Stats of TE model:')
% disp([mean(T_E_3R(range_model)), var(T_E_3R(range_model)),skewness(T_E_3R(range_model)),kurtosis(T_E_3R(range_model))]);
% disp('Stats of hW model:')
% disp([mean(h_W_3R(range_model)), var(h_W_3R(range_model)),skewness(h_W_3R(range_model)),kurtosis(h_W_3R(range_model))]);
% disp('Stats of u model:')
% disp([mean(u_3R(range_model)), var(u_3R(range_model)),skewness(u_3R(range_model)),kurtosis(u_3R(range_model))]);

%% Persistence
range_model2 = 12*20+1:length(I)-12*10;
a3=T_E_obs_new;
cor_te_obs=nan(12,13);
for i=1:12
    for j=0:12
        if i+j<=12
            a=corrcoef(a3(i:12:length(t_obs)),a3(j+i:12:length(t_obs)));
            cor_te_obs(i,j+1)=a(1,2);
        else
            a=corrcoef(a3(i:12:length(t_obs)-12),a3(j+i:12:length(t_obs)));
            cor_te_obs(i,j+1)=a(1,2);
        end
    end
end

% a3=T_E_3R(range_model2);
% cor_te_model=nan(12,13);
% for i=1:12
%     for j=0:12
%         if i+j<=12
%             a=corrcoef(a3(i:12:length(t_obs)),a3(j+i:12:length(t_obs)));
%             cor_te_model(i,j+1)=a(1,2);
%         else
%             a=corrcoef(a3(i:12:length(t_obs)-12),a3(j+i:12:length(t_obs)));
%             cor_te_model(i,j+1)=a(1,2);
%         end
%     end
% end

a3=T_E_3R(range_model2);
cor_te_model=nan(12,13);
for i=1:12
    for j=0:12
        if i+j<=12
            a=corrcoef(a3(i:12:length(range_model2)),a3(j+i:12:length(range_model2)));
            cor_te_model(i,j+1)=a(1,2);
        else
            a=corrcoef(a3(i:12:length(range_model2)-12),a3(j+i:12:length(range_model2)));
            cor_te_model(i,j+1)=a(1,2);
        end
    end
end


a3=T_C_obs_new;
cor_tc_obs=nan(12,13);
for i=1:12
    for j=0:12
        if i+j<=12
            a=corrcoef(a3(i:12:length(t_obs)),a3(j+i:12:length(t_obs)));
            cor_tc_obs(i,j+1)=a(1,2);
        else
            a=corrcoef(a3(i:12:length(t_obs)-12),a3(j+i:12:length(t_obs)));
            cor_tc_obs(i,j+1)=a(1,2);
        end
    end
end
%
% a3=T_C_3R(range_model2);
% cor_tc_model=nan(12,13);
% for i=1:12
%     for j=0:12
%         if i+j<=12
%             a=corrcoef(a3(i:12:length(t_obs)),a3(j+i:12:length(t_obs)));
%             cor_tc_model(i,j+1)=a(1,2);
%         else
%             a=corrcoef(a3(i:12:length(t_obs)-12),a3(j+i:12:length(t_obs)));
%             cor_tc_model(i,j+1)=a(1,2);
%         end
%     end
% end


a3=T_C_3R(range_model2);
cor_tc_model=nan(12,13);
for i=1:12
    for j=0:12
        if i+j<=12
            a=corrcoef(a3(i:12:length(range_model2)),a3(j+i:12:length(range_model2)));
            cor_tc_model(i,j+1)=a(1,2);
        else
            a=corrcoef(a3(i:12:length(range_model2)-12),a3(j+i:12:length(range_model2)));
            cor_tc_model(i,j+1)=a(1,2);
        end
    end
end

figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 15 12])
subplot(2,2,1)
hold on
contourf(0:12,1:12,cor_te_obs,-1:0.1:1,'linestyle','none');
hold on
caxis([-1 1]);
% [cc,hh]=contour(0:12,1:12,cor_te_obs,[0.5 0.5],'linewidth',2,'color','k');
% clabel(cc,hh,'fontsize',9);
set(gca,'xtick',0:1:12,'xticklabel',0:1:12);
set(gca,'ytick',1:12,'yticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});
% xlabel('Lead time (months)','fontsize',9);
ylabel('Initial month','fontsize',9);
set(gca,'fontsize',9);
hold off
colormap(jet);
% colorbar
title('(a) T_E Obs','fontsize',11);

subplot(2,2,2)
hold on
contourf(0:12,1:12,cor_te_model,-1:0.1:1,'linestyle','none');
hold on
caxis([-1 1]);
% [cc,hh]=contour(0:12,1:12,cor_te_model,[0.5 0.5],'linewidth',2,'color','k');
% clabel(cc,hh,'fontsize',9);
set(gca,'xtick',0:1:12,'xticklabel',0:1:12);
set(gca,'ytick',1:12,'yticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});
% xlabel('Lead time (months)','fontsize',9);
% ylabel('Initial month','fontsize',9);
set(gca,'fontsize',9);
hold off
colormap(jet);
% colorbar
title('(b) T_E Model','fontsize',11);

subplot(2,2,3)
hold on
contourf(0:12,1:12,cor_tc_obs,-1:0.1:1,'linestyle','none');
hold on
caxis([-1 1]);
% [cc,hh]=contour(0:12,1:12,cor_tc_obs,[0.5 0.5],'linewidth',2,'color','k');
% clabel(cc,hh,'fontsize',9);
set(gca,'xtick',0:1:12,'xticklabel',0:1:12);
% set(gca,'ytick',1:12,'yticklabel',{'Jan','Feb','Mar','Apr','May','Jun',...
%     'Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',1:12,'yticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});
xlabel('Lead time (months)','fontsize',9);
ylabel('Initial month','fontsize',9);
set(gca,'fontsize',9);
hold off
colormap(jet);
% colorbar
title('(c) T_C Obs','fontsize',11);

subplot(2,2,4)
hold on
contourf(0:12,1:12,cor_tc_model,-1:0.1:1,'linestyle','none');
hold on
caxis([-1 1]);
% [cc,hh]=contour(0:12,1:12,cor_tc_model,[0.5 0.5],'linewidth',2,'color','k');
% clabel(cc,hh,'fontsize',9);
set(gca,'xtick',0:1:12,'xticklabel',0:1:12);
set(gca,'ytick',1:12,'yticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});
xlabel('Lead time (months)','fontsize',9);
% ylabel('Initial month','fontsize',9);
set(gca,'fontsize',9);
hold off
colormap(jet);
title('(d) T_C Model','fontsize',11);
colorbar('eastoutside','position',[.92,0.22,0.01,0.6],'yTick',-1:0.2:1,'fontsize',9,'fontname','arial');

