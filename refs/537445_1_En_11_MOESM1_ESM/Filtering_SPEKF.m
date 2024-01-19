% Chapter 7
% Filtering intermittent time series using the SPEKF 
% The results of the MSM is also provided for comparison
% In generate the true signal, a two-state Markov jump process is used for
% the hidden variable
rng(1) % fix the random number seed
N = 150000; % total number of time stpes
dt = 0.005; % numerical integration time step
dt_obs = 0.5; % observational time step
obs_gap = dt_obs/dt; % the number of numerical integration steps between two observations
ro = 0.5; % another choice is 0.2; % observational noise
% observational operator; note that the one dimensional complex variable u 
% is written as two dimensional real-valued variables when applying the 
% Bayesian update in filtering with real[u] and imag[u]
g = [1,0,0;0,1,0]; 
u = zeros(1,N); % state variable u
gamma = zeros(1,N); % state variable gamma
% model parameters
sigma_u = 0.4;
omega = 2;
f_u = 0; % another choice is f_u = 5;
mu = 0.2;
nu = 0.2;
gamma_plus = 2.27;
gamma_minus = -0.04;
gamma(1) = gamma_plus;
% noise inflation coefficient in the MSM
inflation = 1; % another value is inflation = 3;
% two-state Markov jump process
for i = 2:N
    u(i) = u(i-1) + (-gamma(i-1) + 1i*omega) * u(i-1) * dt + f_u * dt + sigma_u * sqrt(dt) * (randn + 1i*randn)/sqrt(2);
    rd = rand;
    if gamma(i-1) == gamma_minus        
        if rd < mu * dt
            gamma(i) = gamma_plus;
        else
            gamma(i) = gamma_minus;
        end
    else
        if rd < nu * dt
            gamma(i) = gamma_minus;
        else
            gamma(i) = gamma_plus;
        end
    end
end
% calibrating the three parameters in the gamma process of the SPEKF model
% by matching the three statistics with the truth: mean, variance and
% decorrelation time
hat_gamma = (mu * gamma_plus + nu * gamma_minus) / (mu + nu); % mean value of the hidden variable
lag = 1500; % lag of computing the ACF
ACF = autocorr(gamma,lag); % computing the ACF
d_gamma = 1/(-[0:dt:lag*dt]'\log(abs(ACF'))); % the inverse of the decorrelation time is the damping coefficient
sigma_gamma = sqrt(var(gamma) * 2 * d_gamma); % calibrating the noise coefficient
N = 40000; % the total time steps up to which the filtering is carried out
obs_num = length(1:obs_gap:N); % total number of observations given the observational time step
t_obs = [1:obs_gap:N]*dt; % total range of time
u_obs = u(1:obs_gap:N) + ro * (randn(1,obs_num) + 1i * randn(1,obs_num)) / sqrt(2); % observations with noise

u_real = real(u(1:N)); % real part of the observed signal u
u_imag = imag(u(1:N)); % imaginary part of the observed signal u

gamma_truth = gamma(1:N); % the true value of gamma
Num = 2000; % number of ensembles
u0 = zeros(Num,1); % initial ensemble of u
gamma0 = zeros(Num,1); % initial ensemble of gamma
post_mean = zeros(3,obs_num-1); % posterior mean of both u and gamma
post_var = zeros(1,obs_num-1); % posterior variance of gamma
u0_MSM = zeros(Num,1); % initial value of u in the MSM 
post_mean_MSM = zeros(2,obs_num-1); % posterior mean in the MSM
% filtering
% note that analytic formulae can be utilized for forecasting the SPEKF
% model, but here we simply use Monte Carlo simulation in the forecast step
% for the convenience of presenting the results. With a sufficiently large
% number of ensemble size, the sampling error remains negligible
for k = 1:obs_num-1    
    for i = 1:obs_gap % forecast from one observational time step to the next
        u = u0 + (-gamma0 + 1i*omega) .* u0 * dt + f_u * dt + sigma_u * sqrt(dt) * (randn(Num,1) + 1i*randn(Num,1))/sqrt(2);
        gamma = gamma0 - d_gamma * (gamma0 - hat_gamma) * dt + sigma_gamma *sqrt(dt) * randn(Num,1);
        u0 = u;
        gamma0 = gamma;
        
        u_MSM = u0_MSM + (-hat_gamma + 1i*omega) .* u0_MSM * dt + f_u * dt + inflation * sigma_u * sqrt(dt) * (randn(Num,1) + 1i*randn(Num,1))/sqrt(2);
        u0_MSM = u_MSM;
    end
    % SPEKF
    prior_mean = mean([real(u), imag(u), gamma])'; % prior mean 
    prior_cov = cov([real(u), imag(u), gamma]); % prior covariance
    Kalman_Gain = prior_cov * g' * inv(g * prior_cov * g' + ro^2 * eye(2)); % Kalman gain
    post_mean(:,k) = prior_mean + Kalman_Gain * ([real(u_obs(k+1)); imag(u_obs(k+1))] - g * prior_mean); % posterior mean
    post_cov = (eye(3) - Kalman_Gain* g) * prior_cov; % posterior covariance
    post_cov = (post_cov + post_cov')/2; % a numerical trick to guarantee the positive definite property of the covariance matrix
    post_var(k) = post_cov(3,3); % only save the variance in gamma 
    sampling = mvnrnd(post_mean(:,k), post_cov,Num); % generating samplings from posterior states, serving as the initial value for the next assimilation cycle
    u0 = sampling(:,1) + sampling(:,2) * 1i; % initial value of u of the next cycle
    gamma0 = sampling(:,3); % initial value of gamma of the next cycle
    % MSM
    prior_mean_MSM = mean([real(u_MSM), imag(u_MSM)])'; % prior mean
    prior_cov_MSM = cov([real(u_MSM), imag(u_MSM)]); % prior covariance
    Kalman_Gain_MSM = prior_cov_MSM * eye(2) * inv(eye(2) * prior_cov_MSM*eye(2) + ro^2 * eye(2)); % Kalman gain
    post_mean_MSM(:,k) = prior_mean_MSM + Kalman_Gain_MSM * ([real(u_obs(k+1)); imag(u_obs(k+1))] - eye(2) * prior_mean_MSM); % posterior mean
    post_cov_MSM = (eye(2) - Kalman_Gain_MSM* eye(2)) * prior_cov_MSM; % posterior covariance
    post_cov_MSM = (post_cov_MSM + post_cov_MSM') / 2; % a numerical trick to guarantee the positive definite property of the covariance matrix
    sampling = mvnrnd(post_mean_MSM(:,k), post_cov_MSM,Num); % generating samplings from posterior states, serving as the initial value for the next assimilation cycle
    u0_MSM = sampling(:,1) + sampling(:,2) * 1i; % initial value of u of the next cycle        
end
% showing the results
figure
subplot(2,1,1)
hold on
plot(dt:dt:N*dt, u_real,'b','linewidth',2)
plot(t_obs, real(u_obs),'ko','linewidth',2)
plot(t_obs(2:end), post_mean(1,:),'r','linewidth',2)
plot(t_obs(2:end), post_mean_MSM(1,:),'g','linewidth',2)
set(gca,'fontsize',16)
box on
xlim([80,180])
legend('Truth','Obs','SPEKF','MSM')
title('Observed variable Re[u]')
subplot(2,1,2)
hold on
plot(dt:dt:N*dt, gamma_truth,'b','linewidth',2)
plot(t_obs(2:end), post_mean(3,:),'r','linewidth',2)
% showing the uncertainty in the posterior estimate
patch([t_obs(2:end), t_obs(end:-1:2)],[post_mean(3,:) + sqrt(post_var), post_mean(3,end:-1:1) - sqrt(post_var(end:-1:1))],'r','facealpha',0.2,'linestyle','none')
plot([dt,N*dt], [hat_gamma,hat_gamma],'g','linewidth',2)
set(gca,'fontsize',16)
box on
xlim([80,180])
xlabel('t')
title('Unobserved variable \gamma')