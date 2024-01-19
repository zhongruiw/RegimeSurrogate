% Chapter 6
% Computing the approximate response of the SPEKF model
% 1) the qG FDT based on the SPEKF model
% 2) the FDT based on the mean stochastic model (MSM)
% 3) the 'idealized response' but using Monte Carlo simulation rather than 
% the close analytic formula. The Monte Carlo simulation based calculation 
% suffers from the sampling error even for this two dimensional system due 
% to the appearance of the strong non-Gaussian PDF with a one-sided fat tail.  

rng(1) % fix the random number seed
N = 500000; % number of points in time to show the model trajectories
% model parameters
f = 1; 
sigma = 0.5;
sigma_gamma = 1;
d_gamma = 1.3;
gamma_hat = 1;
% forcing perturbation parameters
aa = 1; 
tc = 2;

u = zeros(1,N); % state variable u
gamma = zeros(1,N); % state variable gamma
dt = 0.005; % numerical integration time step
tnum = 2400; % number of points to compute the ACF of u in the SPEKF model, which will be used to calibrate the MSM
% numerical integration of the SPEKF model without adding the perturbation
for i = 2:N     
    u(i) = u(i-1) - gamma(i-1) * u(i-1) * dt + f * dt + sigma * sqrt(dt) * randn; 
    gamma(i) = gamma(i-1) - d_gamma * (gamma(i-1) - gamma_hat) * dt + sigma_gamma * sqrt(dt) * randn; 
end
% showing the time series and statistics of the unperturbed model
figure
subplot(2,4,1:3)
plot(dt:dt:N*dt,u,'b','linewidth',2)
set(gca,'fontsize',12);
box on
title('Time series of u','fontsize',16)
subplot(2,4,4)
[fi,xx]=ksdensity(u);
plot(xx,fi,'b','linewidth',2)
set(gca,'fontsize',12);
box on
title('PDF of u','fontsize',16)
subplot(2,4,[1:3]+4)
plot(dt:dt:N*dt,gamma,'b','linewidth',2)
set(gca,'fontsize',12);
box on
title('Time series of \gamma','fontsize',16)
subplot(2,4,4+4)
[fi,xx]=ksdensity(gamma);
plot(xx,fi,'b','linewidth',2)
set(gca,'fontsize',12);
box on
title('PDF of \gamma','fontsize',16)

umean_att = mean(u); % equilibrium mean of u
uvar_att = var(u); % equilibrium variance of u

ACF = crosscorr(u,u,tnum); % computing the cross correlation
ACF = ACF(tnum+1:end); % computing the ACF
ucorr_att = trapz(0:dt:tnum*dt, ACF); % computing the decorrelation time, which will be used to calibrate the MSM

Sigma = cov([u',gamma']); % covariance of u and gamma
time = 0:dt:tnum*dt; % time window with perturbation
F = 0.1*(tanh(aa * ( time - tc) ) + tanh(aa*tc))/(1+tanh(aa*tc)); % forcing perturbation
vv = Sigma^(-1) * [u - mean(u); gamma - mean(gamma)]; % computing the FDT
vv = vv(1,10000:end); % ignoring the burn-in period
uu = u(10000:end); % ignoring the burn-in period
gg = gamma(10000:end); % ignoring the burn-in period
xcf = crosscorr(vv,uu,tnum); % cross correlation 
xcf = xcf(tnum+1:end);
xcf2 = crosscorr(vv, (uu - mean(uu)).^2, tnum); % cross correlation 
xcf2 = xcf2(tnum+1:end);
resp1 = zeros(1,tnum+1);
resp2 = zeros(1,tnum+1);
for i = 2:tnum+1 % computing the FDT operators
    resp1(i) = trapz(dt:dt:i*dt,xcf(i:-1:1).*F(1:i));
    resp2(i) = trapz(dt:dt:i*dt,xcf2(i:-1:1).*F(1:i));
end
% showing the results
figure
% use the Monte Carlo method to compute the idealized response
time = 0:dt:tnum*dt;
a = 1;
Num = 100000; % ensemble number
u0 = u(10000:4:10000-4+4*Num)'; u0temp = u0; % initial value of u
gamma0 = gamma(10000:4:10000-4+4*Num)'; gamma0temp = gamma0; % initial value of gamma
u_mean = zeros(1,tnum+1); % mean of u
u_var = zeros(1,tnum+1); % variance of u
u_mean(1) = 0; % initial value in mean
u_var(1) = 0; % initial value i variance
for i = 1:tnum+1 % Monte Carlo simulation
    F = 0.1*(tanh(aa * ( i*dt - tc) ) + tanh(aa*tc))/(1+tanh(aa*tc))+f;
    u = u0 + (-gamma0 .* u0 + F) * dt + randn(Num,1) * sqrt(dt) * sigma;
    gamma = gamma0 - d_gamma * (gamma0 - gamma_hat) * dt + sigma_gamma * sqrt(dt) * randn(Num,1); 
    u0 = u;
    gamma0 = gamma;
    u_mean(i) = mean(u);
    u_var(i) = var(u);
end

% calibration the MSM by matching the mean, variance and decorrelation time
% to determine the three parameters
dm = 1/ucorr_att;
sigmam = sqrt(2*dm*uvar_att);
fm = umean_att * dm;

% computing the response using the MSM
% here, Monte Carlo simulation is used, but direct FDT can also be utilized
Num = 100000;
u0 = u0temp;
u_mean2 = zeros(1,tnum+1);
u_var2 = zeros(1,tnum+1);
u_mean2(1) = 0;
u_var2(1) = 0;
for i = 1:tnum+1
    F = 0.1*(tanh(aa * ( i*dt - tc) ) + tanh(aa*tc))/(1+tanh(aa*tc))+fm;
    u = u0 + (-dm .* u0 + F) * dt + randn(Num,1) * sqrt(dt) * sigmam;    
    u0 = u;
    u_mean2(i) = mean(u);
    u_var2(i) = var(u);
end
% showing the results
subplot(3,1,1)
hold on
plot(time, resp1 + u_mean_t0, 'r', 'linewidth',2)
plot(time, u_mean - u_mean(1) + u_mean_t0, 'k', 'linewidth',2)
plot(time, u_mean2 - u_mean2(1) + u_mean_t0, 'g', 'linewidth',2)
set(gca,'fontsize',12)
box on
legend('qG FDT','MC','MSM')
title('(a)  \langle{u}\rangle','fontsize',14,'fontname','times')
subplot(3,1,2)
hold on
plot(time, resp2 + u_var_t0,'r', 'linewidth',2)
plot(time, u_var - u_var(1) + u_var_t0, 'k', 'linewidth',2)
plot(time, u_var2*0 - u_var2(1)*0 + u_var_t0, 'g', 'linewidth',2) % there is no response of variance in the MSM
set(gca,'fontsize',12)
box on
title('(b)  Var(u)','fontsize',14,'fontname','times')
subplot(3,1,3) % time evolution of the forcing; perturbation starts at t = 0.
plot(s,ff,'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(c)  f_u','fontsize',14,'fontname','times')
xlim([0,s(end)])
xlabel('t')
