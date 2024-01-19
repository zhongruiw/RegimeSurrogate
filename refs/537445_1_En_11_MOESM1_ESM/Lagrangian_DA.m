% Chapter 8
% Lagrangian data assimilation with a given number of tracer L
rng(3); % fix the random number seed to reproduce results
sigma_xy = 0.1; % noise in the Lagrangian tracer equations
% L = 20; % comment this out if this code is not run together with
% Lagrangian_DA_Different_L.m
x = zeros(L,N);
y = zeros(L,N);
x(:,1) = rand(L,1)*2*pi;
y(:,1) = rand(L,1)*2*pi;

for i = 2:N % generating the tracer locations
    x(:,i) = x(:,i-1) + exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) * (u_hat(:,i-1) .* transpose(rk(1,:))) * dt + randn(L,1) * sigma_xy * sqrt(dt);
    y(:,i) = y(:,i-1) + exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) * (u_hat(:,i-1) .* transpose(rk(2,:))) * dt + randn(L,1) * sigma_xy * sqrt(dt);
    x(:,i) = mod(x(:,i) + pi, 2*pi) - pi; % periodic boundary conditions
    y(:,i) = mod(y(:,i) + pi, 2*pi) - pi; % periodic boundary conditions
    
end

l = length(k(1,:)); % number of Fourier wavenumbers
if sigma_g ~=0
    sgm = [sigma_g^2 * ones(1,2*l), sigma_B^2 * ones(1,l-1), sigma_x_b^2, sigma_y_b^2];
    dp = [d_g* ones(1,2*l), d_B* ones(1,l-1), d_b, d_b];
    R_eq  = diag(sgm/2*2./dp);
    mu_eq = zeros(3*l+1,1);
    Dim = 3*l+1;
else
    sgm = [sigma_B^2 * ones(1,l-1), sigma_x_b^2, sigma_y_b^2];
    dp = [d_B* ones(1,l-1), d_b, d_b];
    R_eq  = diag(sgm/2*2./dp);
    mu_eq = zeros(l+1,1);
    Dim = l+1;
end
% quantify the uncertainty reduction using relative entropy
Relative_Entropy_Signal = zeros(1,N);
Relative_Entropy_Dispersion = zeros(1,N);
% a matrix used in the filtering formulae
InvBoB = eye(2*L)/sigma_xy/sigma_xy;
mu0 = u_hat(:,1); % initial value of posterior mean
n = length(kk(1,:));
R0 = zeros(n,n); % initial value of posterior covariance
u_post_mean = zeros(n,N); % posterior mean
u_post_mean(:,1) = mu0;
u_post_cov = zeros(n,N); % posterior covariance
u_post_cov(:,1) = diag(R0); % only save the diagonal elements
for i = 2:N
    x0 = x(:,i-1); 
    y0 = y(:,i-1); 
    x1 = x(:,i); 
    y1 = y(:,i); 
    x_diff = x1-x0;
    y_diff = y1-y0;
    % need to take into account the periodic boundary conditions
    x_diff(x_diff > pi)  = x_diff(x_diff>pi)  - 2 * pi;
    x_diff(x_diff < -pi) = x_diff(x_diff<-pi) + 2 * pi;
    y_diff(y_diff > pi)  = y_diff(y_diff>pi)  - 2 * pi;
    y_diff(y_diff < -pi) = y_diff(y_diff<-pi) + 2 * pi;
    % matrix for filtering
    A1 = zeros(2*L,n);
    A1(1:L,:)     = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(L,1) * rk(1,:));
    A1(L+1:2*L,:) = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(L,1) * rk(2,:));
    % update the posterior mean and posterior covariance
    mu = mu0 + (a0 + a1 * mu0) * dt + (R0 * A1') * InvBoB * ([x_diff;y_diff] - A1 * mu0 * dt);
    R = R0 + (a1 * R0 + R0* a1' + Sigma_u*Sigma_u' - (R0*A1') * InvBoB * (R0*A1')')*dt;
    u_post_mean(:,i) = mu;
    u_post_cov(:,i) = diag(R);
    mu0 = mu;
    R0 = R;
    if sigma_g ~=0
        mu_t = mu;
        R_t = R;
    else
        mu_t = mu(2*l+1:end);
        R_t = R(2*l+1:end,2*l+1:end);
    end
    % computing the information gain via relative entropy
    Relative_Entropy_Signal(i) = real(1/2 * ( (mu_t-mu_eq)' / R_eq * (mu_t-mu_eq) ));
    Relative_Entropy_Dispersion(i) = real(1/2 * ( trace(R_t/R_eq) - Dim - log(det( R_t/R_eq)) ));
end
Relative_Entropy_Signal_All = mean(Relative_Entropy_Signal(1000:end));
Relative_Entropy_Dispersion_All = mean(Relative_Entropy_Dispersion(1000:end));

% The following lines are for plotting the results
% figure
% for i = 1:4
%     subplot(2,2,i)
%     if i == 1
%         hold on
%         plot(dt:dt:N*dt, u_hat(1,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, u_post_mean(1,:), 'r', 'linewidth',2)
%         title(['(a) Gravity mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
%     elseif i == 2
%         hold on
%         plot(dt:dt:N*dt, u_hat(Dim_Ug*2+1,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, u_post_mean(Dim_Ug*2+1,:), 'r', 'linewidth',2)
%         title(['(b) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
%     elseif i == 3
%         hold on
%         plot(dt:dt:N*dt, u_hat(6,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, u_post_mean(6,:), 'r', 'linewidth',2)
%         title(['(c) Gravity mode ( ', num2str(kk(1,6)),' , ', num2str(kk(2,6)), ' )'],'fontsize',14)
%     elseif i == 4
%         hold on
%         plot(dt:dt:N*dt, u_hat(Dim_Ug*2+6,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, u_post_mean(Dim_Ug*2+6,:), 'r', 'linewidth',2)
%         title(['(d) GB mode ( ', num2str(kk(1,6)),' , ', num2str(kk(2,6)), ' )'],'fontsize',14)
%     end
%     set(gca,'fontsize',12)
%     box on
%     xlabel('t')
% end
figure
ss = 1;
for i = 2:2:round(T/dt/100)
    
    u = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([0, 2*pi ])
    ylim([0, 2*pi ])
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    box on    
    title(['t = ', num2str(dt*100*(i-1))])    
    hold on
    plot(x(:,1+100*(i-1)),y(:,1+100*(i-1)),'ko')    
    pause(0.1);
    hold off
    ss = ss + 1;
    set(gca,'fontsize',12)
end
figure
for i = 1:9
    subplot(3,3,i)
    hold on
    u = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([0, 2*pi ])
    ylim([0, 2*pi ])
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    box on    
    title(['t = ', num2str(dt*100*(i-1))])
    plot(x(:,1+100*(i-1)),y(:,1+100*(i-1)),'ko','linewidth',2)    
    set(gca,'fontsize',12)
end