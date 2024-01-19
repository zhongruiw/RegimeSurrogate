% Chapter 5.
% Generate the truth signal of the noisy L63 model for the test using ETKF
rng(10) % fix the random number seed
T = 40; % total time length
dt = 0.005; % numerical integration time step 
dt_obs = .5; % observational time step
obs_noise = 3; % observational noise
N = round(T/dt); % total numerical integration steps
N_gap = round(dt_obs/dt); % observational gap: every 'N_gap' numerical integration steps to have one observation
% model parameters
sigma = 10;
rho = 28;
beta = 8/3;
sigma_x = sqrt(2);  
sigma_y = sqrt(12);
sigma_z = sqrt(12);
% the state variables
x_truth = zeros(1,N);
y_truth = zeros(1,N);
z_truth = zeros(1,N);
% initial values
x_truth(1) =  1.5+0.1;
y_truth(1) = -1.5+0.1;
z_truth(1) =  25+0.1;
% numerical integration, generating the true signal
for i = 2:N
    x_truth(i) = x_truth(i-1) + sigma * (y_truth(i-1) - x_truth(i-1)) * dt + sigma_x * sqrt(dt) * randn;
    y_truth(i) = y_truth(i-1) + (x_truth(i-1) * (rho - z_truth(i-1)) - y_truth(i-1)) * dt + sigma_y * sqrt(dt) * randn;
    z_truth(i) = z_truth(i-1) + (x_truth(i-1) * y_truth(i-1) - beta * z_truth(i-1)) * dt + sigma_z * sqrt(dt) * randn;
end
% adding observational noise to the true signal at the observational time
% instants
x_obs = x_truth(1: N_gap: end) + randn(1,N/N_gap) * obs_noise;
y_obs = y_truth(1: N_gap: end) + randn(1,N/N_gap) * obs_noise;
z_obs = z_truth(1: N_gap: end) + randn(1,N/N_gap) * obs_noise;
% showing the results
figure
subplot(3,1,1)
hold on
plot(dt:dt:N*dt, x_truth, 'b', 'linewidth',2);
plot(dt:N_gap*dt:N*dt, x_obs, 'ko', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('x','fontsize',12)
legend('Truth','Obs')
subplot(3,1,2)
hold on
plot(dt:dt:N*dt, y_truth, 'b', 'linewidth',2);
plot(dt:N_gap*dt:N*dt, y_obs, 'ko', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('y','fontsize',12)
subplot(3,1,3)
hold on
plot(dt:dt:N*dt, z_truth, 'b', 'linewidth',2);
plot(dt:N_gap*dt:N*dt, z_obs, 'ko', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('z','fontsize',12)
xlabel('t')