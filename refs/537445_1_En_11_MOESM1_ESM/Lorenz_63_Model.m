% Chapter 1
% Simulate a chaotic system: Lorenz 63 model

T = 40; % total time length
dt = 0.005; % numerical integration time step 
N = round(T/dt); % total numerical integration steps

% model parameters
sigma = 10;
rho = 28;
beta = 8/3;

% the first simulation
% state variables
x1 = zeros(1,N);
y1 = zeros(1,N);
z1 = zeros(1,N);
% initial values
x1(1) =  1.5;
y1(1) = -1.5;
z1(1) =  25;
% model simulation
for i = 2:N
    x1(i) = x1(i-1) + sigma * (y1(i-1) - x1(i-1)) * dt;
    y1(i) = y1(i-1) + (x1(i-1) * (rho - z1(i-1)) - y1(i-1)) * dt;
    z1(i) = z1(i-1) + (x1(i-1) * y1(i-1) - beta * z1(i-1)) * dt;
end

% the second simulation
% state variables
x2 = zeros(1,N);
y2 = zeros(1,N);
z2 = zeros(1,N);
% initial value; differs from the first one by a small amount
epsilon = 0.1;
x2(1) =  1.5 + epsilon;
y2(1) = -1.5 + epsilon;
z2(1) =  25 + epsilon;
% model simulation
for i = 2:N
    x2(i) = x2(i-1) + sigma * (y2(i-1) - x2(i-1)) * dt;
    y2(i) = y2(i-1) + (x2(i-1) * (rho - z2(i-1)) - y2(i-1)) * dt;
    z2(i) = z2(i-1) + (x2(i-1) * y2(i-1) - beta * z2(i-1)) * dt;
end

figure
subplot(3,1,1)
hold on
plot(dt:dt:N*dt, x1, 'b', 'linewidth',2);
plot(dt:dt:N*dt, x2, 'r', 'linewidth',2);
legend('Simulation 1','Simulation 2')
set(gca,'fontsize',12)
box on
title('x','fontsize',16)
subplot(3,1,2)
hold on
plot(dt:dt:N*dt, y1, 'b', 'linewidth',2);
plot(dt:dt:N*dt, y2, 'r', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('y','fontsize',16)
subplot(3,1,3)
hold on
plot(dt:dt:N*dt, z1, 'b', 'linewidth',2);
plot(dt:dt:N*dt, z2, 'r', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('z','fontsize',16)
xlabel('t')


