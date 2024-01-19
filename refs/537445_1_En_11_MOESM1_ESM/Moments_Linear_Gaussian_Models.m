% Chapter 4: 
% Time evolution of the moments for linear Gaussian systems
% dx = (- a x + f) dt + sigma dW  

rng(1); % fix the random number seed
a = 1; % parameter a in the model
f = 1; % parameter f in the model
sigma = 1; % parameter sigma in the model
T = 3; % Total time length
dt = 0.005; % numerical intergration time step
N = round(T/dt); % total points in time
Num = 1000; % the number of samples for the Monte Carlo simulation
x_MC = zeros(Num,N); % point-wise solution from Monte Carlo simulation
mu_MC = zeros(1,N); % associated time evolution of the mean from Monte Carlo simulation
var_MC = zeros(1,N); % associated time evolution of the variance from Monte Carlo simulation
mu_MC(1) = -5; % initial value mean
var_MC(1) = 0.002; % initial value of variance

x_MC(:,1) = randn(Num, 1) * sqrt(var_MC(1)) + mu_MC(1); % initial samples of Monte Carlo simulation based on the given mean and variance
for i = 2:N % time evolution of Monte Carlo
    x_MC(:,i) = x_MC(:,i-1) - a * x_MC(:,i-1) * dt + f * dt + sigma * sqrt(dt) * randn(Num,1); % ensemble solution of SDE
    mu_MC(i) = mean(x_MC(:,i)); % mean 
    var_MC(i) = var(x_MC(:,i)); % variance
end
t = 0: dt: (N-1)*dt; % time points
mu_analytic = mu_MC(1) * exp(- a * t) + f / a * (1 - exp(- a * t) ); % analytic solution of the mean
var_analytic = var_MC(1) * exp(- 2 * a * t) + sigma^2 / 2 / a * (1 - exp(- 2 * a * t)); % analytic solution of the variance
% showing the results
figure
subplot(3, 1, 1) % the ensemble evolution of the solution
plot(t, x_MC, 'g', 'linewidth', 0.5)
set(gca,'fontsize',12)
box on
title('Monte Carlo simulation')
subplot(3, 1, 2) % time evolution of the mean
hold on
plot(t, mu_MC, 'b', 'linewidth', 2)
plot(t, mu_analytic, 'r', 'linewidth', 2)
set(gca,'fontsize',12)
box on
legend('Monte Carlo','Analytic solution')
title('Time evolution of mean')
subplot(3, 1, 3) % time evolution of the variance
hold on
plot(t, var_MC, 'b', 'linewidth', 2)
plot(t, var_analytic, 'r', 'linewidth', 2)
set(gca,'fontsize',12)
box on
title('Time evolution of variance')
xlabel('t')