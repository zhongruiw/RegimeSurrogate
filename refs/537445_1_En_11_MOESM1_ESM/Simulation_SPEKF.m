% Chapter 7
% SPEKF model and its simulation
rng(12) % fix the random number seed
N = 200000; % total numbers of the numerical integration time steps
dt = 0.005; % numerical integration time step
u = zeros(1,N); % observed state variable u
gamma = zeros(1,N); % hidden variable: stochastic damping gamma
omega = zeros(1,N); % hidden variable: stochasitc phase omega
b = zeros(1,N); % hidden variable: stochastic forcing b

% model parameters
sigma_u = 1;

d_gamma = 1;
gamma_hat = 1;
sigma_gamma = 1.2;

d_omega = 0.5;
omega_hat = 5;
sigma_omega = 1.5;

d_b = 0.25;
b_hat = 0;
sigma_b = 15;
% simulate the model
for i = 2:N
    u(i) = u(i-1) + ( (- gamma(i-1) + 1i * omega(i-1)) * u(i-1) + b(i-1)) * dt + sigma_u * sqrt(dt) * (randn + 1i*randn)/sqrt(2);
    gamma(i) = gamma(i-1) - d_gamma * (gamma(i-1) - gamma_hat) * dt + sigma_gamma * sqrt(dt) * randn;
    omega(i) = omega(i-1) - d_omega * (omega(i-1) - omega_hat) * dt + sigma_omega * sqrt(dt) * randn;
    b(i) = b(i-1) - d_b * (b(i-1) - b_hat) * dt + sigma_b * sqrt(dt) * (randn + 1i * randn)/sqrt(2);
end
% showing the anti-damping phases
gamma_index = find(gamma<=0);
x = dt:dt:N*dt;
gamma_negative = gamma(gamma_index);
x_negative = x(gamma_index);
gamma_index2 = gamma_index(1);
j = 2;
for i = 2:length(gamma_index)
    if gamma_index(i) - gamma_index(i-1)>1.5
        gamma_index2 = [gamma_index2, gamma_index(i-1), gamma_index(i)];
        j = j + 1;
    end
end
gamma_index2 = reshape(gamma_index2(1:end-1),2,[]);

% showing the values of the phase omega outside one standard deviation from the mean
k1 = 2;
omega_index = find(omega>omega_hat + sqrt(sigma_omega^2/2/d_omega));
omega_index2 = omega_index(1);
for i = 2:length(omega_index)
    if omega_index(i) - omega_index(i-1)>1.5
        omega_index2 = [omega_index2, omega_index(i-1), omega_index(i)];
        k1 = k1 + 1;
    end
end
omega_index2 = reshape(omega_index2(1:end-1),2,[]);

k2 = 2;
omega_index = find( omega<omega_hat - sqrt(sigma_omega^2/2/d_omega));
omega_index3 = omega_index(1);
for i = 2:length(omega_index)
    if omega_index(i) - omega_index(i-1)>1.5
        omega_index3 = [omega_index3, omega_index(i-1), omega_index(i)];
        k2 = k2 + 1;
    end
end
omega_index3 = reshape(omega_index3(1:end-1),2,[]);

% showing the values of forcing f outside one standard deviation from the mean
k3 = 2;
b_index = find(real(b)>b_hat + sqrt(sigma_b^2/2/d_b));
b_index2 = b_index(1);
for i = 2:length(b_index)
    if b_index(i) - b_index(i-1)>1.5
        b_index2 = [b_index2, b_index(i-1), b_index(i)];
        k3 = k3 + 1;
    end
end
b_index2 = reshape(b_index2(1:end-1),2,[]);

k4 = 2;
b_index = find(real(b)<b_hat - sqrt(sigma_b^2/2/d_b));
b_index3 = b_index(1);
for i = 2:length(b_index)
    if b_index(i) - b_index(i-1)>1.5
        b_index3 = [b_index3, b_index(i-1), b_index(i)];
        k4 = k4 + 1;
    end
end
b_index3 = reshape(b_index3(1:end-1),2,[]);

% showing the results
figure
subplot(4,1,1)
hold on
plot(dt:10*dt:N*dt, real(u(1:10:end)),'b','linewidth',2)
plot([dt,N*dt],[0,0],'k--')
box on
set(gca,'fontsize',16)
xlim([180,310])
title('Time series of SPEKF model')
subplot(4,1,2)
hold on
plot(dt:10*dt:N*dt, gamma(1:10:end),'b','linewidth',2)
for jj = 1:j-2
    plot([x(gamma_index2(1,jj):gamma_index2(2,jj))],[gamma(gamma_index2(1,jj):gamma_index2(2,jj))],'r','linewidth',2)
end
plot([dt,N*dt],[0,0],'k--')
box on
set(gca,'fontsize',16)
xlim([180,310])
subplot(4,1,3)
hold on
plot(dt:10*dt:N*dt, omega(1:10:end),'b','linewidth',2)
for kk = 1:k1-2
    plot([x(omega_index2(1,kk):omega_index2(2,kk))],[omega(omega_index2(1,kk):omega_index2(2,kk))],'g','linewidth',2)
end
for kk = 1:k2-2
    plot([x(omega_index3(1,kk):omega_index3(2,kk))],[omega(omega_index3(1,kk):omega_index3(2,kk))],'c','linewidth',2)
end
box on
set(gca,'fontsize',16)
xlim([180,310])
subplot(4,1,4)
hold on
plot(dt:10*dt:N*dt, real(b(1:10:end)),'b','linewidth',2)
box on
set(gca,'fontsize',16)
xlim([180,310])
for kk = 1:k3-2
    plot([x(b_index2(1,kk):b_index2(2,kk))],[real(b(b_index2(1,kk):b_index2(2,kk)))],'g','linewidth',2)
end
for kk = 1:k4-2
    plot([x(b_index3(1,kk):b_index3(2,kk))],[real(b(b_index3(1,kk):b_index3(2,kk)))],'c','linewidth',2)
end
xlabel('t')