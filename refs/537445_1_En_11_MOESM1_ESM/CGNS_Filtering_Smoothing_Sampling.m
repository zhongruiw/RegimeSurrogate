% Chapter 8
% Using CGNS for filtering, smoothing and conditional sampling
% The test model here is a nonlinear dyad model, where
% u is the observed variable and v is unobserved and needs to be recovered


rng(123) % fix the random number seed for repeating experiments
% Parameters
N = 100000; % total number of numerical integration time steps
dt = 0.005; % numerical integration time step;  
 
u_truth = zeros(1,N); % state variable u
v_truth = zeros(1,N); % state variable v

% Dimension of the hidden variables
% the code here is written for general cases where the dimension of the
% hidden variables is not necessarily to be one
Dim = 1; % Dimension of the hidden variables
Dim2 = 1; %Dim2 is Dim^2, the dimension of the covariance matrix;

% Model parameters
F_u = 1;
F_v = 1;
d_v = 0.8;
d_u = 0.8;
c = 1.2;
sigma_v = 2;
sigma_u = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generating the truth signal %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the model forward to generate the true signals
for i = 2:N
    u_truth(i) = u_truth(i-1) + ( ( - d_u + c * v_truth(i-1)) * u_truth(i-1) + F_u) * dt + sqrt(dt) * randn * sigma_u;
    v_truth(i) = v_truth(i-1) + (- d_v * v_truth(i-1) - c * u_truth(i-1)^2 + F_v) * dt + sqrt(dt) * randn * sigma_v;
end
 
% Plot the true signal and statistics
figure
subplot(2,4,1:3) % plot trajectory
plot(dt:dt:N*dt,u_truth,'b','linewidth',1)
box on
set(gca,'fontsize',16) 
title('(a) A trajectory of the observed variable u','fontsize',14)
subplot(2,4,4)
[fi,xx] = ksdensity(u_truth(1:10:end));
fi_G = normpdf(xx,mean(u_truth(1:10:end)),std(u_truth(1:10:end)));
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'k--','linewidth',2)
legend('Truth','Gaussian fit')
box on
set(gca,'fontsize',16) 
title('(b) PDF of u','fontsize',14)
subplot(2,4,5:7) % plot trajectory
plot(dt:dt:N*dt,v_truth,'b','linewidth',2)
box on
set(gca,'fontsize',16) 
title('(c) A trajectory of the hidden variable v','fontsize',14)
xlabel('t')
subplot(2,4,8)
[fi,xx] = ksdensity(v_truth(1:10:end));
fi_G = normpdf(xx,mean(v_truth(1:10:end)),std(v_truth(1:10:end)));
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'k--','linewidth',2)
legend('Truth','Gaussian fit')
box on
set(gca,'fontsize',16) 
title('(d) PDF of v','fontsize',14)
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save the results of the posterior mean and posterior covariance
gamma_mean_trace = zeros(Dim,N);
gamma_cov_trace = zeros(Dim2,N);
% initial value of the posterior mean and posterior covariance
gamma_mean0 = zeros(Dim,1);
gamma_cov0 = eye(Dim)*0.01; % The initial covariance is arbitary for filtering but it must be nonzero to avoid singularity in smoothing
gamma_cov_trace(:,1) = reshape(gamma_cov0,Dim2,1); % save the initial covariance

for i = 2:N
    % matrices and vectors in the conditional Gaussian framework
    u0 = u_truth(i-1);
    u = u_truth(i);
    a1 = - d_v;
    a0 = - c * u0^2 + F_v;
    b1 = sigma_v;
    A0 = -d_u * u0 + F_u;
    A1 = c * u0;
    invBoB =  1/sigma_u/sigma_u;

    % time evolution of the posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (u-u0 - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
    
    % save the posterior statistics
    gamma_mean_trace(:,i) = gamma_mean;
    gamma_cov_trace(:,i) = gamma_cov;

    % update
    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Smoothing and backward sampling %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_s = zeros(Dim,N); % Save the posterior mean in smoothing
R_s = zeros(Dim2,N); % Save the posterior covariance in smoothing
C = zeros(Dim2,N); % an intermediate matrix used in smoothing
 
% Smoothing is backward
% Intial values for smoothing (at the last time instant)
Num = 3; % the number of sampled trajectories 
rng(2); % fix the random number seed to reproduce the results
mu_s(:,end) = gamma_mean; % save the initial value of the smoother mean 
R_s(:,end) = reshape(gamma_cov,Dim2,1); % save the initial value of the smoother covariance
R_s_temp = gamma_cov; % an intermediate matrix
Y_Sampling = zeros(Num,N); % Save for the backward sampled trajectory
rd_Y = randn(Num,N); % pre-generated random numbers
for i = N-1:-1:1
    % Matrices and vectors in the conditional Gaussian smoothing and
    % backward sampling
    u0 = u_truth(:,i);
    a1 = - d_v;
    a0 = - c * u0^2 + F_v;    
    gamma_cov = reshape(gamma_cov_trace(:,i),Dim,Dim); % filter covariance is needed as the input of smoothing formula 
    C_temp = gamma_cov * (eye(Dim) + a1 * dt)' * (b1 * b1' * dt + (eye(Dim) + a1 * dt) * gamma_cov * (eye(Dim) + a1 * dt)')^(-1);
    C(:,i) = reshape(C_temp,Dim2,1);
    mu_s(:,i) = gamma_mean_trace(:,i) + C_temp * (mu_s(:,i+1) - a0 * dt - ( eye(Dim) + a1 * dt) * gamma_mean_trace(:,i)); % update the smoother mean
    R_s_temp = reshape(R_s(:,i+1),Dim,Dim);
    R_s_temp = gamma_cov + C_temp * (R_s_temp - (eye(Dim) + a1 * dt) * gamma_cov * (eye(Dim) + a1 * dt)' - b1 * b1 * dt) * C_temp';   
    R_s(:,i) = reshape(R_s_temp,Dim2,1); % update the smoother covariance (this line and the above two lines)
    Y_Sampling(:,i) = Y_Sampling(:,i+1) + (-a0 - a1 * Y_Sampling(:,i+1)) * dt + b1 * b1' * inv(gamma_cov) * (gamma_mean_trace(:,i) ... 
        - Y_Sampling(:,i+1)) * dt + b1 * rd_Y(:,i) * sqrt(dt); % Backward sampling; the sampled trajectory has random noise
end





figure
subplot(4,4,1:3) % plot trajectory
h0=plot(dt:dt:N*dt,u_truth,'g','linewidth',1);
box on
set(gca,'fontsize',16) 
title('(a) Time series','fontsize',16)
grid on
grid(gca,'minor')
xlim([275,290])
set(gca,'xtick',275:290)
subplot(4,4,4)
[fi,xx] = ksdensity(u_truth(1:10:end));
hold on
plot(xx,fi,'g','linewidth',2)
box on
set(gca,'fontsize',16) 
title('(b) PDFs','fontsize',16)
subplot(4,4,5:7) % plot trajectory
hold on
h1=plot(dt:dt:N*dt,v_truth,'b','linewidth',2);
h2=plot(dt:dt:N*dt,gamma_mean_trace,'r','linewidth',2);
gamma_upper = gamma_mean_trace+2*sqrt(gamma_cov_trace);
gamma_lower = gamma_mean_trace-2*sqrt(gamma_cov_trace);
tt = dt:dt:N*dt;
patch([tt,tt(end:-1:1)],[gamma_lower,gamma_upper(end:-1:1)],'r','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',16) 
xlim([275,290])
set(gca,'xtick',275:290)
grid on
grid(gca,'minor')
subplot(4,4,8)
[fi,xx] = ksdensity(v_truth(1:10:end));
[fi2,xx] = ksdensity(gamma_mean_trace(1:10:end),xx);
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi2,'r','linewidth',2)
box on
set(gca,'fontsize',16) 
subplot(4,4,9:11) % plot trajectory
hold on
plot(dt:dt:N*dt,v_truth,'b','linewidth',2)
h3=plot(dt:dt:N*dt,mu_s,'c','linewidth',2);
gamma_upper = mu_s+2*sqrt(R_s);
gamma_lower = mu_s-2*sqrt(R_s);
tt = dt:dt:N*dt;
patch([tt,tt(end:-1:1)],[gamma_lower,gamma_upper(end:-1:1)],'c','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',16) 
xlim([275,290])
set(gca,'xtick',275:290)
grid on
grid(gca,'minor')
subplot(4,4,12)
[fi,xx] = ksdensity(v_truth(1:10:end));
[fi3,xx] = ksdensity(mu_s(1:10:end),xx);
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi3,'c','linewidth',2)
box on
set(gca,'fontsize',16) 
subplot(4,4,13:15) % plot trajectory
hold on
plot(dt:dt:N*dt,v_truth,'b','linewidth',2)
h4=plot(dt:dt:N*dt,Y_Sampling(1,:),'color',[0.1,0.1,0.1],'linewidth',2);
h5=plot(dt:dt:N*dt,Y_Sampling(2,:),'color',[0.4,0.4,0.4],'linewidth',2);
h6=plot(dt:dt:N*dt,Y_Sampling(3,:),'color',[0.7,0.7,0.7],'linewidth',2);
legend([h0,h1,h2,h3,h4,h5,h6],'Truth of u','Truth of v','Filter mean','Smoother mean','Sampling 1','Sampling 2','Sampling 3')
box on
set(gca,'fontsize',16) 
xlabel('t')
xlim([275,290])
set(gca,'xtick',275:290)
grid on
grid(gca,'minor')
subplot(4,4,16)
[fi,xx] = ksdensity(v_truth(1:10:end));
[fi4,xx] = ksdensity(Y_Sampling(1,1:10:end),xx);
[fi5,xx] = ksdensity(Y_Sampling(2,1:10:end),xx);
[fi6,xx] = ksdensity(Y_Sampling(3,1:10:end),xx);
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi4,'color',[0.1,0.1,0.1],'linewidth',2)
plot(xx,fi5,'color',[0.4,0.4,0.4],'linewidth',2)
plot(xx,fi6,'color',[0.7,0.7,0.7],'linewidth',2)
box on
set(gca,'fontsize',16) 