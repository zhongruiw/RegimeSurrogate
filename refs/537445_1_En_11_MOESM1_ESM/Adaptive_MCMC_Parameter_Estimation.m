% Chapter 9
% Adaptive MCMC for the cubic model
% du = (f + au + bu^2 - cu^3) dt + sigma dW
% The reference of the adaptive MCMC is Vihola (2012). 

rng(10); % fix the random number seed

T = 500; % total length of the observed time series
dt = 5e-3; % numerical integration time step
N = round(T/dt); % total number of time steps within [0,T]
u_truth = zeros(1,N); % true signal (observations)


Lag = 2000; % lag in computing the ACF, nothing to do with MCMC, just for showing the model property
% true parameters (choose one of the two regimes)
% a_truth = -2; b_truth = 2.5; c_truth = 1; f_truth = 0.0; sigma_truth = 1; %highly skewed regime
a_truth = -2; b_truth = 0; c_truth = 0; f_truth = 0.5; sigma_truth = 0.5; % Gaussian regime 


% generate the true signal
for i = 2:N
    u_truth(i) = u_truth(i-1) + (f_truth + a_truth * u_truth(i-1) + b_truth * u_truth(i-1)^2 - c_truth * u_truth(i-1)^3) * dt + sqrt(dt) * randn * sigma_truth;
end
u_obs = u_truth;
% plotting the true signal
figure
subplot(3,5,1:3)
hold on
plot(dt:dt:N*dt, u_truth,'b', 'linewidth',2)
set(gca,'fontsize',12)
box on
title('Observed signal from true system')
subplot(3,5,4)
ACF_truth = autocorr(u_truth, Lag);
hold on
plot(0:dt:Lag*dt,ACF_truth,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF')
subplot(3,5,5)
[fi, xx] = ksdensity(u_truth); fi_truth_u = fi; xx_truth_u = xx;
fi_G = normpdf(xx,mean(u_truth),sqrt(var(u_truth)));
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'k--','linewidth',2)
set(gca,'fontsize',12)
box on
legend('Truth','Gaussian fit')
title('PDF')
 

%%% Estimate the parameters using adaptive MCMC

% initial guess of the parameters
a = 2; b = -2; c = 2; f = 2; sigma = 2;

% MCMC setup
KK = 10000; % total MCMC iteraction steps
KK_BurnIn = 500; % burn-in period
FLAG = zeros(1,KK); % FLAG(k) = 1 or 0, representing accept or reject the proposed value at the k-th step
LIKELI = zeros(1,KK); % store the likelihood at each iteraction
alpha_star = 0.25; % target acceptance ratio
num_param = 5; % total number of parameters
S = .1 * eye(num_param); % initial guess of S for adaptive MCMC
alpha_k = 1; % acceptence ratio
% store all the trace plot [f, a, b, c, log(sigma)]; this will guarantee
% the noise coefficient sigma > 0
Parameters_all = zeros(num_param,KK); % store the trace plot
Parameters_all(:,1) = [f, a, b, c, sigma];
eta0 = 25; % a tuning parameter for the decaying of eta used in the adaptive MCMC


Likelihood1 = -500000;% assigning the initial likelihood, no need to be accurate

% MCMC iteraction
for kk = 2:KK 
    if mod(kk,1000) == 0
        disp(['    Iter step kk = ', num2str(kk)])
        disp('    f,         a,         b,         c,       \sigma')
        disp([Parameters_all(1,kk-1), Parameters_all(2,kk-1), Parameters_all(3,kk-1), Parameters_all(4,kk-1), exp(Parameters_all(5,kk-1))])
    end
    
    if kk >= KK_BurnIn % applying the adaptive MCMC after the burn-in period
        eta = min(1, num_param * eta0 * kk^(-2/3));
        S = chol( S * (eye(num_param) + eta * (exp(alpha_k) - alpha_star)*(U*U')/(U'*U) )*S');
        S = S';
    end
    
    if kk < KK_BurnIn-1 % no adaptive for the burn-in period
        Parameters_all(:,kk) = Parameters_all(:,kk-1) + 0.05*randn(num_param,1);
    else % applying the adaptive MCMC after the burn-in period
        U = randn(num_param,1) * (1 + norm(Parameters_all(:,kk-1))^2)^((-num_param+1)/2);
        Parameters_all(:,kk) = Parameters_all(:,kk-1) + S * U;
    end
    
    
    % save the results 
    f = Parameters_all(1,kk);
    a = Parameters_all(2,kk);
    b = Parameters_all(3,kk);
    c = Parameters_all(4,kk);
    sigma = exp(Parameters_all(5,kk));

    % computing the log-likelihood using a local linear Gaussian
    % approximation
    effective_damping = a + 2*b*u_obs(1:end-1) - 3*c*u_obs(1:end-1).^2;
    u_mean = u_obs(1:end-1) .* exp( effective_damping * dt) + ( 1- exp(effective_damping * dt)) * f ./ (-effective_damping);
    u_var = (1 - exp(2 * effective_damping * dt)) / 2 ./ (-effective_damping) * sigma^2;

    likelitemp = log( normpdf(u_obs(2:end),u_mean,sqrt(u_var)) ); % likelihood based on the Gaussian approximation
    likelitemp(isnan(likelitemp)) = 0; likelitemp(isinf(likelitemp)) = -1000; % if the likelihood is -infity, then set it to be a large negative number
    Likelihood2 = sum( likelitemp );


    arate = (Likelihood2 - Likelihood1); % computing the acceptance ratio

    alpha_k = min(0, arate); % accept or reject
    if alpha_k > log(rand)
        flag = 1; % accept
        Likelihood1 =  Likelihood2;  
    else
        flag = 0; % reject
        Parameters_all(:,kk) = Parameters_all(:,kk-1); % using the parameter values in the previous step
    end
    FLAG(kk) = flag;
    LIKELI(kk) = Likelihood1;


end


f_all = Parameters_all(1,:);
a_all = Parameters_all(2,:);
b_all = Parameters_all(3,:);
c_all = Parameters_all(4,:);
Parameters_all(5,:) = exp(Parameters_all(5,:));
sigma_all = Parameters_all(5,:);

disp('acceptance ratio: expected,     actual')
disp([alpha_star,    sum(FLAG(KK_BurnIn:end))/(KK-KK_BurnIn)])



% trace plot and posterior pdf

for i = 1:num_param
    if i == 1
        param_truth = f_truth;
    elseif i == 2
        param_truth = a_truth;
    elseif i == 3
        param_truth = b_truth;
    elseif i == 4
        param_truth = c_truth;
    elseif i == 5
        param_truth = sigma_truth;
    end
    subplot(3,num_param,num_param+i)
    hold on
    plot(1:KK, Parameters_all(i,:),'b','linewidth',2)
    plot([1,KK],[param_truth,param_truth],'k--','linewidth',2)
    box on
    set(gca,'fontsize',12)
    xlim([1,KK])
    if i == 1
        title('f')
    elseif i == 2
        title('a')
    elseif i == 3
        title('b')
    elseif i == 4
        title('c')
    elseif i == 5
        title('\sigma')
    end
    xlabel('k')
end
% simulting the model using the estimated parameters
% the estimated parameters are defined by taking the average after the
% burn-in period
u_est = zeros(1,N); 
f_est = mean(Parameters_all(1,KK-3*KK_BurnIn:end));
a_est = mean(Parameters_all(2,KK-3*KK_BurnIn:end));
b_est = mean(Parameters_all(3,KK-3*KK_BurnIn:end));
c_est = mean(Parameters_all(4,KK-3*KK_BurnIn:end));
sigma_est = mean(Parameters_all(5,KK-3*KK_BurnIn:end));
for i = 2:N
    u_est(i) = u_est(i-1) + (f_est + a_est * u_est(i-1) + b_est * u_est(i-1)^2 - c_est * u_est(i-1)^3) * dt + sqrt(dt) * randn * sigma_est;
end


subplot(3,5,[1:3]+10)
hold on
plot(dt:dt:N*dt, u_est,'r', 'linewidth',2)
set(gca,'fontsize',12)
box on
title('Trajectory of u using the estimated parameters')
xlabel('t')
subplot(3,5,4+10)
ACF_est = autocorr(u_est, Lag);
hold on
plot(0:dt:Lag*dt,ACF_est,'r','linewidth',2)
plot(0:dt:Lag*dt,ACF_truth,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF')
xlabel('t')
subplot(3,5,5+10)
[fi, xx] = ksdensity(u_est); fi_est_u = fi; xx_est_u = xx;
fi_G = normpdf(xx,mean(u_est),sqrt(var(u_est)));
hold on
plot(xx,fi,'r','linewidth',2)
% plot(xx,fi_G,'k--','linewidth',2)
plot(xx_truth_u,fi_truth_u,'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('PDF')
legend('Model with estimation','Truth')
xlabel('u')
