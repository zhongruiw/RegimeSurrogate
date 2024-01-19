% Chapter 7
% The three dimensional ENSO model
% Note that the code only shows the model simulation but not the
% observational data. The observational data can be downloaded from the
% resources provided in the main text
rng(10) % fix the random number seed
figure
N = 248000; % total number of time steps in a long simulation
dt = 0.005; % numerical integration time step
Point_40year = round(40/dt); % number of time steps corresponding to 40 years, which is the length of the observations
Total_Periods = round(N/Point_40year); % total segments, each having 40 years
Te = zeros(1,N); % state variable: SST in the eastern Pacific
Hw = zeros(1,N); % state variable: thermocline depth in the western Pacific
tau = zeros(1,N); % state variable: wind bursts
% model parameters
d_T = 1.5; 
d_H = 1.5;
d_tau = 4;
omega = -1.5;
alpha_T = 1;
alpha_H = -0.4;
sigma_T = 0.8;
sigma_H = 0.8;
Lag = round(5/dt); % lag used to compute the ACF
for i = 2:N % numerical intergration
    sigma_tau = 4.5 * (tanh(Te(i-1))+1)+4;
    Te(i) = Te(i-1) + (-d_T * Te(i-1) - omega * Hw(i-1) + alpha_T * tau(i-1)) * dt + sigma_T * randn * sqrt(dt);
    Hw(i) = Hw(i-1) + (-d_H * Hw(i-1) + omega * Te(i-1) + alpha_H * tau(i-1)) * dt + sigma_H * randn * sqrt(dt);
    tau(i) = tau(i-1) + (-d_tau * tau(i-1)) * dt + sigma_tau * randn * sqrt(dt);
end
% showing the results of the model
for i = 1:3
    if i == 1
        variable = Te;
    elseif i == 2
        variable = Hw;
    else
        variable = tau;
    end
    subplot(3,6,[1:4]+(i-1)*6)
    plot([dt:dt:N*dt]+1800, variable,'r','linewidth',2) % time series
    box on
    set(gca,'fontsize',16);
    if i == 1
        ylabel('^oC')
    elseif i == 2
        ylabel('15m')
    else
        ylabel('m/s');        
    end
    xlabel('Year');
    if i == 1
        title('(a) Trajectories')
    end
    xlim([1982,2020])
    subplot(3,6,5+(i-1)*6)
    ACF = autocorr(variable, Lag);
    plot(0:dt:Lag*dt, ACF,'r','linewidth',2) % ACF
    box on
    set(gca,'fontsize',16)
    if i == 1
        title('(b) ACFs')
    end
    xlabel('Year')
    ACF_all = zeros(Total_Periods - 1,Lag+1);
    for j = 1:Total_Periods-1
        ACF_all(j,:) = autocorr(variable([1:Point_40year] + j * Point_40year), Lag);
    end
    ACF_mean = mean(ACF_all);
    ACF_var = var(ACF_all);
    ACF_upper = ACF_mean + 2*sqrt(ACF_var);
    ACF_lower = ACF_mean - 2*sqrt(ACF_var); 
    hold on % uncertainty of the ACF 
    patch([0:dt:Lag*dt,Lag*dt:-dt:0],[ACF_lower,ACF_upper(end:-1:1)],'r','facealpha',0.15,'linestyle','none')
    subplot(3,6,6+(i-1)*6)
    [pdf,xx] = ksdensity(variable); % PDF
    plot(xx,pdf,'r','linewidth',2);
    box on
    set(gca,'fontsize',16);
    if i == 1
        title('(c) PDFs')
        xlabel('^oC')
    elseif i == 2
        xlabel('15m')
    else
        xlabel('m/s')
    end
    pdf_all = zeros(Total_Periods - 1,100);
    for j = 1:Total_Periods-1
        [pdf_all(j,:),xx] = ksdensity(variable([1:Point_40year] + j * Point_40year), xx);
    end
    pdf_mean = mean(pdf_all);
    pdf_var = var(pdf_all);
    pdf_upper = pdf_mean + 2*sqrt(pdf_var);
    pdf_lower = pdf_mean - 2*sqrt(pdf_var);pdf_lower(pdf_lower<0)=0;
    hold on % uncertainty of the PDF 
    patch([xx,xx(end:-1:1)],[pdf_lower,pdf_upper(end:-1:1)],'r','facealpha',0.15,'linestyle','none')
end