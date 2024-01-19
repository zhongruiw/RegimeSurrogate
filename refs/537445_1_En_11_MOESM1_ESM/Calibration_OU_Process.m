% Chapter 7
% Calibration a complex OU process
rng(1); % fix the random number seed
N = 100000; % total number of time steps
dt = 0.005; % numerical integration time step
u = zeros(1,N); % state variable 
% true parameters in the complex OU process
gamma = 1; % real-valued damping
f = 1 + 1*1i; % complex-valued forcing 
sigma = 1; % real-valued noise
omega = -3; % real-valued phase
for i = 2:N % numerical integration
    u(i) = u(i-1) + (-gamma + 1i * omega) * u(i-1) * dt + f * dt + sigma/sqrt(2) * sqrt(dt) * (randn + 1i * randn);
end

figure
% showing the observed time series of u
subplot(2,1,1)
hold on
plot(dt:dt:N*dt, real(u), 'b', 'linewidth',2)
plot(dt:dt:N*dt, imag(u), 'r', 'linewidth',2)
set(gca,'fontsize',12)
box on
title('Observed true signal')
legend('Real[u]','Imag[u]')
xlabel('t')

% model calibration
Lag = 1000; % lag for computing the ACF
r = autocorr(real(u),Lag); % ACF of the real part of u
r2 = crosscorr(real(u),imag(u),Lag); % cross-correlation between the real and imaginary parts of u
tt = 0:dt:Lag*dt; % time interval to plot the ACF or cross-correlation function
F = @(x,data)exp(-x(1)*data).*sin(x(2)*data); % the ansatz of the cross-correlation 
% fitting the cross-correlation function
x0 = [0.5,0.5];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,tt,r2(Lag+1:end)); % numerical fit
% showing the ACF and the cross-correlation
subplot(2,1,2)
hold on
plot(0:dt:Lag*dt, r, 'b', 'linewidth', 2)
plot(0:dt:Lag*dt, r2(Lag+1:end), 'r', 'linewidth', 2)
box on
set(gca,'fontsize', 12)
title('ACF and Cross-correlation function')
legend('ACF','Cross-correlation function')
xlabel('t')

gamma_est = x(1); % estimated damping gamma
omega_est = x(2); % estimated phase omega

m = mean(u); % mean of u
E = var(u); % variance of u
%%%% Note: gamma^2 + omega^2 = T^2 + theta^2;
T = gamma_est/(gamma_est^2+omega_est^2); % parameter T, the real part of the decorrelation time
theta = omega/(gamma_est^2+omega_est^2); % parameter theta, the imaginary part of the decorrelation time
f_est = m * (T-1i*theta)/(T^2+theta^2); % estimated forcing f
sigma_est = sqrt(2*E*T/(T^2+theta^2)); % estimated noise coefficient sigma
disp('Truth parameters:')
disp('\gamma,     \omega,     f,     \sigma')
disp([gamma,omega,f,sigma]);
disp('Estimated values:')
disp('\gamma,     \omega,     f,     \sigma')
disp([gamma_est,omega_est,f_est,sigma_est]);