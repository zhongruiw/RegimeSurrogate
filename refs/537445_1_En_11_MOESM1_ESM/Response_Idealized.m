% Chapter 6
% Computing the idealized response of the SPEKF model
% Here the idealized response is computed by solving the time evolution of 
% the statistics via closed analytic formulae of the SPEKF model. The
% details of the closed analytic formulae can be found in Gershgorin et.
% al., 2010, JCP.

N = 5000; % total number of time integration steps
dt = 0.005; % integration time step
% parameters in the SPEKF model
d_gamma = 1.3;
gamma_hat = 1;
sigma_u = 0.8;
sigma_gamma = 1;
nn = 50; % save the results at every 'nn' steps
u_mean = zeros(1,N/nn); % time evolution of the mean
u_var = zeros(1,N/nn); % time evolution of the variance
% auxiliary variables in solving the time evolution of the statistics in
% the SPEKF model 
A_temp1 = zeros(1,N/nn);
B_temp1 = zeros(1,N/nn);
C_temp1 = zeros(1,N/nn);
AB_temp1 = zeros(1,N/nn);
ff = zeros(1,N/nn); % forcing perturbation
u0 = 2; u_mean(1) = u0; % initial value of u
gamma0 = gamma_hat; % initial value of gamma
% parameters in the ramp-type perturbation
aa = 1;
tc = 2;

ii = 1; % an index to save the data

% intergration via closed analytic formula
for i = 2:nn:N
    if mod(i,500) == 2
        disp(i)
    end
    s = 0:dt:(i-1)*dt;
    t = s(end);
    if i <=2500
        F = 1;       
    else
        aa = 1;
        tc = 2;
        F = [ones(1,2500),0.1*(tanh(aa * ( s(2501:end)-dt*2500 - tc) ) + tanh(aa*tc))/(1+tanh(aa*tc))+1];
    end
    Jst = 1/d_gamma * ( exp(-d_gamma * s) - exp(-d_gamma * t)) * (gamma0 - gamma_hat);
    VarJst = -sigma_gamma^2/d_gamma^3 * ( 1 + d_gamma * (s-t) + exp(-d_gamma * (s+t)) .* (-1 - exp(2 * d_gamma * s) + cosh(d_gamma * (s-t))) );
    u_mean(ii) = exp(-gamma_hat * t) * u0 * exp(- Jst(1) + 1/2* VarJst(1)) + trapz(s, exp(-gamma_hat * (t-s)) .* F .* exp(- Jst + 1/2* VarJst));
    A = exp(-2 * gamma_hat * t) * u0^2 * exp(- 2 * Jst(1) + 2 * VarJst(1));
    
    ss = ones(i,1) * s; rr = ss';
    Jst_temp_all = 1/d_gamma * ( exp(-d_gamma * ss) - exp(-d_gamma * t)) * (gamma0 - gamma_hat);
    Jrt_temp_all = 1/d_gamma * ( exp(-d_gamma * rr) - exp(-d_gamma * t)) * (gamma0 - gamma_hat);
    VarJst_temp_all = -sigma_gamma^2/d_gamma^3 * ( 1 + d_gamma * (ss-t) + exp(-d_gamma * (ss+t)) .* (-1 - exp(2 * d_gamma * ss) + cosh(d_gamma * (ss-t))) );
    VarJrt_temp_all = -sigma_gamma^2/d_gamma^3 * ( 1 + d_gamma * (rr-t) + exp(-d_gamma * (rr+t)) .* (-1 - exp(2 * d_gamma * rr) + cosh(d_gamma * (rr-t))) );
    CovJstJrt_temp_all = - sigma_gamma^2/2/d_gamma^3 * ( exp(-d_gamma * (t-ss)) - exp(-d_gamma * (t-rr)) + exp(-d_gamma * (t+ss)) - exp(-d_gamma * (t+rr)) ...
        - 1 + exp(-d_gamma * (ss-rr)) - exp(-2*d_gamma * ss) + exp(-d_gamma * (ss + rr)) ) + VarJst_temp_all;
    for j = 1:i
        for k = j+1:i
            CovJstJrt_temp_all(k,j) = CovJstJrt_temp_all(j,k);
        end
    end
    
    
    B_temp_all = exp(-gamma_hat * (2*t - ss - rr)) .* exp(-Jst_temp_all - Jrt_temp_all + 1/2 * VarJst_temp_all + 1/2 * VarJrt_temp_all + CovJstJrt_temp_all);
    B_temp = B_temp_all; 
    
 
    B = trapz(s, trapz(s,  B_temp .* (F'*F) ));
 
    Jst_all = 1/d_gamma * ( exp(-d_gamma * s) - exp(-d_gamma * t)) * (gamma0 - gamma_hat);
    Jt0t_all = 1/d_gamma * ( ones(1,i) - exp(-d_gamma * t)) * (gamma0 - gamma_hat);
    VarJst_all  = -sigma_gamma^2/d_gamma^3 * ( 1 + d_gamma * (s-t) + exp(-d_gamma * (s+t)) .* (-1 - exp(2 * d_gamma * s) + cosh(d_gamma * (s-t))) );
    VarJt0t_all = -sigma_gamma^2/d_gamma^3 * ( 1 + d_gamma * (-t) + exp(-d_gamma * (t)) .* (-2 + cosh(d_gamma * (-t))) ) * ones(1,i);
    covJst_all = VarJst_all - sigma_gamma^2/2/d_gamma^3 * ( exp(-d_gamma*(t-s)) - exp(-d_gamma*(t)) + exp(-d_gamma*(t+s)) - exp(-d_gamma*(t)) - 1 + exp(-d_gamma*(s)) - exp(-2*d_gamma*(s)) + exp(-d_gamma*(s))  );


    C = sigma_u^2 * trapz(s, exp(-2 * gamma_hat * (t-s)) .* exp(- 2 * Jst + 2 * VarJst)); 
    AB_temp = - Jt0t_all - Jst_all + 1/2 * VarJt0t_all + 1/2 * VarJst_all + covJst_all;
    AB = exp(-gamma_hat * t) * trapz(s, exp(-gamma_hat*(t-s)) .* F .* u0 .* exp(AB_temp)  );
    u_var(ii) = A + B + C + 2*AB -u_mean(ii)^2;
    A_temp1(ii) = A;
    B_temp1(ii) = B;
    C_temp1(ii) = C;
    AB_temp1(ii) = AB;
    ii = ii + 1;
end



 
% put the forcing and its perturbation into a vector
ii = 1;
for i = 2:nn:N
    if mod(i,500) == 2
        disp(i)
    end
    s = 0:dt:(i-1)*dt;
    t = s(end);
    if i <=2500
        ff(ii) = 1;       
    else
        aa = 1;
        tc = 2;
        ff(ii) = 0.1*(tanh(aa * ( (i-2500)*dt - tc) ) + tanh(aa*tc))/(1+tanh(aa*tc))+1;
    end
    ii = ii + 1;
end


s = [2*dt:nn*dt:(N-1)*dt]-12.5; % rescale the time to make 0 as the starting point of the perturbation

% showing the results
figure
subplot(3,1,1) % time evolution of the mean
plot(s,u_mean,'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(a)  \langle{u}\rangle','fontsize',14,'fontname','times')
xlim([s(1),s(end)])
subplot(3,1,2) % time evolution of the variance
plot(s,u_var,'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(b)  Var(u)','fontsize',14,'fontname','times')
xlim([s(1),s(end)])
subplot(3,1,3) % time evolution of the forcing; perturbation starts at t = 0.
plot(s,ff,'b','linewidth',2)
set(gca,'fontsize',12)
box on
title('(c)  f_u','fontsize',14,'fontname','times')
xlim([s(1),s(end)])
xlabel('t')
