% Chapter 3
% Euler-Maruyama and Milstein schemes for solving the geometric Brownian 
% motion dx = mu X dt + sigma X dW

% parameters
mu = -0.5;
sigma = 2.5;
N = 1000; % total number of integration points
dt = 0.001; % numerical integration time step
X_EM = zeros(1,N); X_EM(1) = 1; % solution of utilizing the Euler-Maruyama scheme
X_M = zeros(1,N); X_M(1) = 1; % solution of utilizing the Milstein scheme
for i = 2:N
    rd = randn; % random number
    X_EM(i) = X_EM(i-1) + mu * X_EM(i-1) * dt + sigma * X_EM(i-1) * sqrt(dt) * rd; % Euler-Maruyama
    X_M(i)  = X_M(i-1)  + mu * X_M(i-1)  * dt + sigma * X_M(i-1)  * sqrt(dt) * rd ...
        + 1/2 * sigma^2 * X_M(i-1) * (rd^2 * dt - dt); % Milstein
end
figure 
hold on
plot(dt:dt:N*dt, X_EM, 'b', 'linewidth',2)
plot(dt:dt:N*dt, X_M,  'r', 'linewidth',2)
legend('Euler-Maruyama','Milstein')
box on
set(gca,'fontsize',12)
title('Model trajectories', 'fontsize',16)