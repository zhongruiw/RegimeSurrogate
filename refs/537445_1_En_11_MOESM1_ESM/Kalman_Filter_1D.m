% Chapter 5
% Kalman filter in one dimension

rng(1); % fix the random number seed
figure
N = 100000; % total number of time steps
dt = 0.005; % numerical integration time step
u = zeros(1,N); % state variable 
% parameters in the linear Gaussian SDE
gamma = 0.4; 
omega_0 = 1;
f_0 = 2;
f_1 = 0;
omega_1 = 0;
sigma = 1;
g = 1; % observational operator
for k = 1:3 % considering three different cases with different observational time step (dt_obs) and the observations noise (r_o) 
    if k == 1
        dt_obs = 0.5;
        r_o = 0.5;
    elseif k == 2
        dt_obs = 0.5;
        r_o = 3;
    elseif k == 3
        dt_obs = 3;
        r_o = 0.5;
    end
    gap = round(dt_obs/dt); % number of numerical integration time steps between every two observations
    for i = 2:N % generating the true signal: integrating the model forward
        u(i) = u(i-1) + ( (-gamma + 1i * omega_0) * u(i-1) + f_0 + f_1 * exp(1i * omega_1 * dt * i)) * dt + sigma * sqrt(dt) * (randn + 1i*randn)/sqrt(2);
    end
    u_truth = u(1:gap:N); % true values at the observational time instants
    L = length(u_truth);
    u_obs = u_truth + (randn(1,L) + 1i * randn(1,L)) / sqrt(2) * sqrt(r_o); % polluted observations
    % introducing variables to save the forecast (pred) and analysis (filter) results 
    u_pred_mean = zeros(1,L); 
    u_pred_var = zeros(1,L);
    u_filter_mean = zeros(1,L);
    u_filter_var = zeros(1,L);
    Kalman_Gain = zeros(1,L); % Kalman gain
    % the for loop repeats the run from one observational time instant to
    % the next one 
    for i = 2:L
        % computing the forecast solution using the analytic expressions
        u_pred_mean(i) = u_filter_mean(i-1) * exp( (-gamma + 1i * omega_0) * dt_obs ) + f_0/(gamma-1i*omega_0) * (1 - exp( (-gamma + 1i * omega_0) * dt_obs) ) ...
            + f_1 * exp(1i*omega_1 * (i*dt_obs+dt))/(gamma + 1i * (-omega_0 + omega_1)) * (1 - exp( -(gamma + 1i * omega_1 - 1i * omega_0) * dt_obs)); 
        u_pred_var(i) = u_filter_var(i-1) * exp(-2 * gamma * dt_obs) + sigma^2/2/gamma * (1 - exp(-2 * gamma * dt_obs));
        % computing the Kalman gain
        Kalman_Gain(i) = g * u_pred_var(i)/(r_o + g^2 * u_pred_var(i));
        % computing the posterior solution
        u_filter_mean(i) = (1-Kalman_Gain(i)*g) * u_pred_mean(i) + Kalman_Gain(i) * u_obs(i);
        u_filter_var(i) = (1-Kalman_Gain(i)*g) * u_pred_var(i);
    end
    tt = dt:dt_obs:N*dt;
    % computing the root-mean-squared error (RMSE) for the forecast and
    % filtering results based on the solutions at the observational time
    % instants, comparing with the true values there
    RMSE_pred = sqrt(sum((u_truth - u_pred_mean).* conj((u_truth - u_pred_mean)))/L);
    RMSE_filter = sqrt(sum((u_truth - u_filter_mean).* conj((u_truth - u_filter_mean)))/L);
    % showing the results
    subplot(3,1,k)
    hold on
    plot(dt:dt:N*dt,u,'b','linewidth',2)
    plot(dt:dt_obs:N*dt,u_obs,'ko','linewidth',2)
    plot(dt:dt_obs:N*dt,u_pred_mean,'g','linewidth',2)
    plot(dt:dt_obs:N*dt,u_filter_mean,'r','linewidth',2)   
    % showing the 95% confidence interval using the prior/posterior standard deviation
    patch([tt,tt(end:-1:1)], [u_pred_mean+2*sqrt(u_pred_var),u_pred_mean(end:-1:1)-2*sqrt(u_pred_var(end:-1:1))],'g','facealpha',0.2,'linestyle','none');
    patch([tt,tt(end:-1:1)], [u_filter_mean+2*sqrt(u_filter_var),u_filter_mean(end:-1:1)-2*sqrt(u_filter_var(end:-1:1))],'r','facealpha',0.2,'linestyle','none');    
    xlim([230,280])
    box on
    set(gca,'fontsize',16)
    if k == 1
        legend('Truth','Obs','Forecast mean','Filter mean')
    end
    if k == 3
        xlabel('t')
    end
    ylim([-4,4])
    text(250,-3.3,['Kalman gain = ', num2str(Kalman_Gain(end)),';  RMSE forecast = ',num2str(RMSE_pred), ';  RMSE filter = ',num2str(RMSE_filter)],'fontsize',16)
    if k == 1
        title(['Case (a): \Delta{t}^{obs} = ',num2str(dt_obs),'; r^o = ',num2str(r_o)],'fontsize',16)
    elseif k == 2
        title(['Case (b): \Delta{t}^{obs} = ',num2str(dt_obs),'; r^o = ',num2str(r_o)],'fontsize',16)
    else
        title(['Case (c): \Delta{t}^{obs} = ',num2str(dt_obs),'; r^o = ',num2str(r_o)],'fontsize',16)
    end
end