% Chapter 9
% Parameter estimation using DA via direct method
% Test model: noisy Lorenz 63
% Using the CGNS as the testbed

rng(3) % fix the random number seed 
T = 30; % total length of observations in time
dt = 0.005; % numerical integration time step
N = round(T/dt); % total number of numerical integration time steps
rho = 28; sgm = 10; beta = 8/3; % parameters
x = zeros(1,N);x(1) = 1; % state variable x
y = zeros(1,N); % state variable y
z = zeros(1,N); % state variable z
sgm_x = 15; % noise coefficient sigma_x
sgm_y = 15; % noise coefficient sigma_y
sgm_z = 15; % noise coefficient sigma_z
for i = 2:N % numerical integration of the true system to generate the observational time series
    x(i) = x(i-1) + sgm * ( y(i-1) - x(i-1) ) * dt + sqrt(dt) * sgm_x * randn;
    y(i) = y(i-1) + ( x(i-1) * ( rho - z(i-1) ) - y(i-1) ) * dt + sqrt(dt) * sgm_x * randn;
    z(i) = z(i-1) + ( x(i-1) * y(i-1) - beta * z(i-1) ) * dt + sqrt(dt) * sgm_z * randn;
end



gamma_mean_trace = zeros(3,N);% posterior mean
gamma_cov_trace = zeros(6,N);% posterior covariance; the posterior covariance is 3x3 but only 6 entries are independent due to the symmetry


% the matrices and vectors in the Y = (parameters) equations
% they are all zeros in the direct method
Sigma_II = zeros(3,3);
a0 = zeros(3,1);
a1 = zeros(3,3);

gamma_mean0 = zeros(3,1); % initial value of the posterior mean
gamma_cov0 = 2 * eye(3); % initial value of the posterior covariance
gamma_cov_trace(1:3,1) = 2; % save the initial covariance

invBoB = inv(diag([sgm_x^2,sgm_y^2,sgm_z^2])); % inverse of the noise matrix square, used for the CGNS filtering formula

% using data assimilation for parameter estimation
for i = 2:N

    % observations
    u = [x(i); y(i); z(i)];
    u0 = [x(i-1); y(i-1); z(i-1)];
    % matrix and vector in the observed equations (i.e., L63 model)
    A0 = [0; - x(i-1) * z(i-1) - y(i-1); x(i-1) * y(i-1)];
    A1 = diag([y(i-1) - x(i-1), x(i-1), -z(i-1)]);


    % updating the posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (u-u0 - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + Sigma_II * Sigma_II' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     

    % save the results
    gamma_mean_trace(:,i) = gamma_mean;
    gamma_cov_trace(1,i) = gamma_cov(1,1);
    gamma_cov_trace(2,i) = gamma_cov(2,2);
    gamma_cov_trace(3,i) = gamma_cov(3,3);
    gamma_cov_trace(4,i) = gamma_cov(1,2);
    gamma_cov_trace(5,i) = gamma_cov(1,3);
    gamma_cov_trace(6,i) = gamma_cov(2,3);


    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;

end
% showing the results
figure
subplot(3,3,1)
plot(dt:10*dt:T,x(1:10:end),'b','linewidth',2)
set(gca,'fontsize',12)
ylim([-30,30]);xlim([0,T]);
box on
title('(a) Sample trajectory of x','fontsize',12)
subplot(3,3,4)
plot(dt:10*dt:T,y(1:10:end),'b','linewidth',2)
set(gca,'fontsize',12)
ylim([-30,30]);xlim([0,T]);
box on
title('(b) Sample trajectory of y','fontsize',12)
subplot(3,3,7)
plot(dt:10*dt:T,z(1:10:end),'b','linewidth',2)
set(gca,'fontsize',12)
ylim([0,60]);xlim([0,T]);
box on
title('(c) Sample trajectory of z','fontsize',12)
xlabel('t')

subplot(3,3,2)
hold on
plot(dt:10*dt:T,gamma_mean_trace(1,1:10:end),'b','linewidth',2)
plot([dt,T],[sgm,sgm],'--k','linewidth',2)
patch([dt:10*dt:T,T:-10*dt:dt],[gamma_mean_trace(1,1:10:end)+2*sqrt(gamma_cov_trace(1,1:10:end)), gamma_mean_trace(1,end:-10:1)-2*sqrt(gamma_cov_trace(1,end:-10:1))],'b','facealpha',0.3,'linestyle','none')
set(gca,'fontsize',12)
box on
title('(d) Estimation of \sigma','fontsize',12);xlim([0,T]);
subplot(3,3,5)
hold on
plot(dt:10*dt:T,gamma_mean_trace(2,1:10:end),'b','linewidth',2)
plot([dt,T],[rho,rho],'--k','linewidth',2)
patch([dt:10*dt:T,T:-10*dt:dt],[gamma_mean_trace(2,1:10:end)+2*sqrt(gamma_cov_trace(2,1:10:end)), gamma_mean_trace(2,end:-10:1)-2*sqrt(gamma_cov_trace(2,end:-10:1))],'b','facealpha',0.3,'linestyle','none')
set(gca,'fontsize',12)
box on
title('(e) Estimation of \rho','fontsize',12);xlim([0,T]);
subplot(3,3,8)
hold on
plot(dt:10*dt:T,gamma_mean_trace(3,1:10:end),'b','linewidth',2)
plot([dt,T],[beta,beta],'--k','linewidth',2)
patch([dt:10*dt:T,T:-10*dt:dt],[gamma_mean_trace(3,1:10:end)+2*sqrt(gamma_cov_trace(3,1:10:end)), gamma_mean_trace(3,end:-10:1)-2*sqrt(gamma_cov_trace(3,end:-10:1))],'b','facealpha',0.3,'linestyle','none')
set(gca,'fontsize',12)
box on
title('(f) Estimation of \beta','fontsize',12);xlim([0,T]);
legend('Estimation','Truth')
xlabel('t')

subplot(3,3,3)
plot(dt:10*dt:T,gamma_cov_trace(1,1:10:end),'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(g) Uncertainty in estimating \sigma','fontsize',12);xlim([0,T]);
subplot(3,3,6)
plot(dt:10*dt:T,gamma_cov_trace(2,1:10:end),'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(h) Uncertainty in estimating \rho','fontsize',12);xlim([0,T]);
subplot(3,3,9)
plot(dt:10*dt:T,gamma_cov_trace(3,1:10:end),'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(i) Uncertainty in estimating \beta','fontsize',12);xlim([0,T]);
xlabel('t')