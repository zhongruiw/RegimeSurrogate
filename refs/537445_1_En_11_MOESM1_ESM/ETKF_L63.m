% Chapter 5.
% ETKF for the noisy Lorenz 63 model
% True signal is generated from the L63_TrueSignal.m, which should run
% first
rng(12) % fix the random number generator
% Dim_obs can take values 1, 2 or 3.
Dim_obs = 1; % Dimension of observations. 1 means x; 2 means x y; 3 means x y z
Dim_system = 3; % Dimension of the L63 model, which is always 3.
Ens_Num = 10; % ensemble size
ini_cov = 2; % initial uncertainty
Ro = obs_noise^2 * eye(Dim_obs);  %observational noise covariance matrix
r = 0; % noise inflation coefficient
% initial ensembles 
x_sample = x_truth(1) + sqrt(ini_cov) * randn(1, Ens_Num);
y_sample = y_truth(1) + sqrt(ini_cov) * randn(1, Ens_Num);
z_sample = z_truth(1) + sqrt(ini_cov) * randn(1, Ens_Num);

posterior_mean_save = zeros(Dim_system, N/N_gap); % posterior mean
posterior_var_save = zeros(Dim_system, N/N_gap); % posterior covariance (only the diagonal entries)
posterior_mean_save(:,1) = [x_truth(1); y_truth(1); z_truth(1)];
v_all = [x_obs;y_obs;z_obs];
v_all = v_all(1:Dim_obs,:); % put all the observations into one matrix  
% observational operator
if Dim_obs == 1 
    G = zeros(1,3); G(1,1) = 1; % observing only x
elseif Dim_obs == 2
    G = zeros(2,3); G(1,1) = 1; G(2,2) = 1; % observing x and y
else
    G = eye(3); % observing x, y and z
end


% sequential ETKF  
for ij = 2:N/N_gap
    % integrating to the next observational time instant
    for i = 2:N_gap
        x_sample_new = x_sample + sigma * (y_sample - x_sample) * dt + sigma_x * sqrt(dt) * randn(1,Ens_Num);
        y_sample_new = y_sample + (x_sample .* (rho - z_sample) - y_sample) * dt + sigma_y * sqrt(dt) * randn(1,Ens_Num);
        z_sample_new = z_sample + (x_sample .* y_sample - beta * z_sample) * dt + sigma_z * sqrt(dt) * randn(1,Ens_Num);
        x_sample = x_sample_new;
        y_sample = y_sample_new;
        z_sample = z_sample_new;
    end
    u_prior = [x_sample; y_sample; z_sample]; % collecting the ensembles
    u_mean_prior = mean(u_prior,2); % ensemble prior mean
    U = u_prior - u_mean_prior * ones(1, Ens_Num); % residual of the ensembles around the mean
    V = G * U; % residual of the observations around the mean; the dimension depends on the number of the observational variables
    J = (Ens_Num - 1) / (1 + r) * eye(Ens_Num) + V' / Ro * V;    
    J = (J + J')/2; % to eliminate the round off numerical error and guarantee J is symmetric
    [X, Gamma] = eig(J); % eigen decomposition
    x = J \ V' / Ro * ( v_all(:,ij) - G * mean(u_prior,2) ); % ETKF
    u_mean_posterior = u_mean_prior + U * x; % posterior mean
    T = sqrt(Ens_Num-1) * X * Gamma^(-1/2) * X';% transform matrix
    U_perturb_posterior = U * T; % posterior perturbation matrix
    u_posterior = u_mean_posterior * ones(1, Ens_Num) + U_perturb_posterior; % posterior ensembles
    posterior_mean_save(:,ij) = u_mean_posterior; % save the posterior mean estimates
    x_sample = u_posterior(1,:);
    y_sample = u_posterior(2,:);
    z_sample = u_posterior(3,:); 
end

% plotting the truth, obs, filter mean
figure
for i = 1:3
    subplot(3,1,i)
    if i == 1
        variable_truth = x_truth;
    elseif i == 2
        variable_truth = y_truth;
    elseif i == 3
        variable_truth = z_truth;
    end
    hold on
    plot(dt:dt:N*dt, variable_truth, 'b', 'linewidth',2);
    if i == 1
        plot(dt:N_gap*dt:N*dt, x_obs, 'ko', 'linewidth',2);
    elseif i == 2
        if Dim_obs >= 2
            plot(dt:N_gap*dt:N*dt, y_obs, 'ko', 'linewidth',2);
        end
    elseif i == 3
        if Dim_obs == 3
            plot(dt:N_gap*dt:N*dt, z_obs, 'ko', 'linewidth',2);
        end
    end
    plot(dt:N_gap*dt:N*dt, posterior_mean_save(i,:), 'm', 'linewidth',2);
    if i == 1
        legend('Truth','Obs','Filter mean (ETKF)')
    end
    box on
    set(gca,'fontsize',12)
end

% computing the skill scores: RMSE and Corr
Truth_all = [x_truth(1:N_gap:end); y_truth(1:N_gap:end); z_truth(1:N_gap:end);];
RMSE = zeros(1,Dim_system);
Corr = zeros(1,Dim_system);
for i = 1:Dim_system
    v1 = Truth_all(i,:);
    v2 = posterior_mean_save(i,:);
    RMSE(i) = sqrt(sum((v1-v2).^2)/length(v1))/std(v1);
    Corr_temp = corrcoef(v1,v2);
    Corr(i) = Corr_temp(1,2);
end
disp('ETKF results')
disp(['RMSE:  ', num2str(RMSE)])
disp(['Corr:  ', num2str(Corr)])