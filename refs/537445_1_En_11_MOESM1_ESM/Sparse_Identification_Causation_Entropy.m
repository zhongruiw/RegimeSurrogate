% Chapter 9
% Sparse identification of a complex dynamical system via information 
% theory (causation entropy)
% The test model used here is the noisy Lorenz 63 model

rng(10) % fix the random number seed to reproduce results
T = 100; % total length of the observational time
dt = 0.001; % numerical integration time step 
N = round(T/dt); % total number of the numerical integration steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generating the true signal %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% true model parameters
sigma = 10;
rho = 28;
beta = 8/3;

sigma_x = 1;  
sigma_y = 1;
sigma_z = 1;

% state variables of the true signal
x_truth = zeros(1,N);
y_truth = zeros(1,N);
z_truth = zeros(1,N);
x_truth(1) = 0;
y_truth(1) = 0;
z_truth(1) = 0;

% numerical integration to generate the true time series
for i = 2:N
    x_truth(i) = x_truth(i-1) + sigma * (y_truth(i-1) - x_truth(i-1)) * dt + sigma_x * sqrt(dt) * randn;
    y_truth(i) = y_truth(i-1) + (x_truth(i-1) * (rho - z_truth(i-1)) - y_truth(i-1)) * dt + sigma_y * sqrt(dt) * randn;
    z_truth(i) = z_truth(i-1) + (x_truth(i-1) * y_truth(i-1) - beta * z_truth(i-1)) * dt + sigma_z * sqrt(dt) * randn;
end
Lag = 1000; % lag used to compute the ACF
% showing the true signal (serving as the observations) and the associated
% statistics (ACF and PDF)
figure
subplot(3,5,1:3)
hold on
plot(dt:dt:N*dt, x_truth, 'b', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('x (true signal)','fontsize',12)
subplot(3,5,4)
ACF = autocorr(x_truth,Lag);
plot(0:dt:Lag*dt,ACF,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of x')
subplot(3,5,5)
[fi,xx] = ksdensity(x_truth);
plot(xx,fi,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of x')
subplot(3,5,[1:3]+5)
hold on
plot(dt:dt:N*dt, y_truth, 'b', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('y','fontsize',12)
subplot(3,5,4+5)
ACF = autocorr(y_truth,Lag);
plot(0:dt:Lag*dt,ACF,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of y')
subplot(3,5,5+5)
[fi,xx] = ksdensity(y_truth);
plot(xx,fi,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of y')
subplot(3,5,[1:3]+10)
hold on
plot(dt:dt:N*dt, z_truth, 'b', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('z','fontsize',12)
xlabel('t')
subplot(3,5,4+10)
ACF = autocorr(z_truth,Lag);
plot(0:dt:Lag*dt,ACF,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of z')
subplot(3,5,5+10)
[fi,xx] = ksdensity(z_truth);
plot(xx,fi,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Sparse identification via Causation entropy %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first thing is to build a library of candidate functions
% Here, all the linear and quadratic nonlinear functions are included in
% the library. More functions can be included if needed.
% These functions are: x, y, z, xy, yz, zx, x^2, y^2, z^2
L = 9; % the total number of candidate functions
% Put all candidate functions into a vector
All_Candidates = [x_truth; 
y_truth; 
z_truth;
x_truth .* y_truth;
y_truth .* z_truth;
z_truth .* x_truth;
x_truth.^2;
y_truth.^2;
z_truth.^2];
% true model structure corresponding to the arrangement of the candidate functions
% x, y, z, xy, yz, zx, x^2, y^2, z^2
True_Model_Structure = ...
[1  1  0  0  0  0  0  0  0
 1  1  0  0  0  1  0  0  0
 0  0  1  1  0  0  0  0  0];

% delete the last point on the time series as the forward Euler scheme is
% utilized in the numerical calculation 
All_Candidates = All_Candidates(:,1:end-1); 
% computing the time derivative (the left hand side of the equations)
x_derivative = (x_truth(2:end)-x_truth(1:end-1))/dt;
y_derivative = (y_truth(2:end)-y_truth(1:end-1))/dt;
z_derivative = (z_truth(2:end)-z_truth(1:end-1))/dt;
% CEM: causation entropy matrix
% three rows correspond to x, y and z
% L columns correspond to the candidate functions
CEM = zeros(3,L);
% computing the CEM using the information theory based on the Gaussian 
% approximation
for k = 1:3
    if k == 1
        All_Variables = [All_Candidates; x_derivative];
    elseif k == 2
        All_Variables = [All_Candidates; y_derivative];
    elseif k == 3
        All_Variables = [All_Candidates; z_derivative];
    end
    All_Cov = cov(All_Variables');
    % computing the causation entropy for each candidate function
    for i = 1:L
        RXY = All_Cov([1:i-1,i+1:end],[1:i-1,i+1:end]);
        RY = All_Cov([1:i-1,i+1:end-1],[1:i-1,i+1:end-1]);
        RXYZ = All_Cov;
        RYZ = All_Cov(1:L,1:L);
        % Gaussian approximation
        CEM(k,i) = 1/2 * ( log(det(RXY)) - log(det(RY)) - log(det(RXYZ)) + log(det(RYZ)) );
    end
end
CEM_Original = CEM;
threshold = 0.0003; % set up a threshold
CEM_indicator = (CEM_Original>threshold); % a 0/1 matrix as an indicator for the model structure
disp('Candidate functions:')
disp('x, y, z, xy, yz, zx, x^2, y^2, z^2')
disp('Model structure (truth):')
disp(num2str(True_Model_Structure))
disp('Model structure (identified model):')
disp(num2str(CEM_indicator))
% below is parameter estimation (after seeing the CEM indicator)
% the parameter estimation is based on a maximum likelihood estimator via a
% quadratic optimization method
% H and g are for including the physics constraints 
H = zeros(1,7); H(5) = 1; H(7) = 1; g = 0;
Theta = zeros(7,1); % parameters in the deterministic part
Sigma = zeros(3,3); % noise coefficients
for ss = 1:10
    for i = 1:N-1
        M = zeros(3,7);
        M(1,1) = All_Candidates(1,i);
        M(1,2) = All_Candidates(2,i);
        M(2,3) = All_Candidates(1,i);
        M(2,4) = All_Candidates(2,i);
        M(2,5) = All_Candidates(6,i);
        M(3,6) = All_Candidates(3,i);
        M(3,7) = All_Candidates(4,i);
        M = M * dt;
        z = [x_truth(i+1);y_truth(i+1);z_truth(i+1)];
        s = [x_truth(i);y_truth(i);z_truth(i)];
        Sigma = Sigma + (z - M * Theta - s) * (z - M * Theta - s)';
    end
    Sigma = Sigma/(N-1).*eye(3);
    invSigma = inv(Sigma);
    D = zeros(7,7);
    c = zeros(7,1);
    for i = 1:N-1
        M = zeros(3,7);
        M(1,1) = All_Candidates(1,i);
        M(1,2) = All_Candidates(2,i);
        M(2,3) = All_Candidates(1,i);
        M(2,4) = All_Candidates(2,i);
        M(2,5) = All_Candidates(6,i);
        M(3,6) = All_Candidates(3,i);
        M(3,7) = All_Candidates(4,i);
        M = M * dt;
        z = [x_truth(i+1);y_truth(i+1);z_truth(i+1)];
        s = [x_truth(i);y_truth(i);z_truth(i)];
        D = D + M' * invSigma * M;
        c = c + M' * invSigma * (z - s);
    end
    D = D/(N-1);
    c = c/(N-1);
    Theta1 = D\c;
    lambda = (H/D*H')^(-1)*(H/D*c-g);
    Theta1_Constraint = D\(c-H'*lambda);
end
Sigma1 = Sigma;


x1 = zeros(1,N);
y1 = zeros(1,N);
z1 = zeros(1,N);
x1(1) = 0;
y1(1) = 0;
z1(1) = 0;
rng(10)
for i = 2:N
    x1(i) = x1(i-1) + (Theta1_Constraint(1) * x1(i-1) + Theta1_Constraint(2) * y1(i-1)) * dt + sqrt(Sigma1(1,1)/dt) * sqrt(dt) * randn;
    y1(i) = y1(i-1) + (Theta1_Constraint(3) * x1(i-1) + Theta1_Constraint(4) * y1(i-1) + Theta1_Constraint(5) * z1(i-1) * x1(i-1)) * dt + sqrt(Sigma1(2,2)/dt) * sqrt(dt) * randn;
    z1(i) = z1(i-1) + (Theta1_Constraint(6) * z1(i-1) + Theta1_Constraint(7) * x1(i-1) * y1(i-1)) * dt + sqrt(Sigma1(3,3)/dt) * sqrt(dt) * randn;
end
figure
subplot(3,5,1:3)
hold on
plot(dt:dt:N*dt, x1, 'r', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('x (identified model with a low threshold)','fontsize',12)
subplot(3,5,4)
ACF = autocorr(x1,Lag);
plot(0:dt:Lag*dt,ACF,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of x')
subplot(3,5,5)
[fi,xx] = ksdensity(x1);
plot(xx,fi,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of x')
subplot(3,5,[1:3]+5)
hold on
plot(dt:dt:N*dt, y1, 'r', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('y','fontsize',12)
subplot(3,5,4+5)
ACF = autocorr(y1,Lag);
plot(0:dt:Lag*dt,ACF,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of y')
subplot(3,5,5+5)
[fi,xx] = ksdensity(y1);
plot(xx,fi,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of y')
subplot(3,5,[1:3]+10)
hold on
plot(dt:dt:N*dt, z1, 'r', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('z','fontsize',12)
xlabel('t')
subplot(3,5,4+10)
ACF = autocorr(z1,Lag);
plot(0:dt:Lag*dt,ACF,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of z')
subplot(3,5,5+10)
[fi,xx] = ksdensity(z1);
plot(xx,fi,'r','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of z')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is the result using a different threshold. So the model structure
% is also different
% Previously, there were 7 parameter, now one parameter is further
% eliminated
threshold = 0.0005;
CEM_indicator = (CEM_Original>threshold);
H = zeros(1,6); H(4) = 1; H(6) = 1; g = 0;
M = zeros(3,6);
Theta = zeros(6,1);
Sigma = zeros(3,3);
for ss = 1:1
    for i = 1:N-1
        M = zeros(3,6);
        M(1,1) = All_Candidates(1,i);
        M(1,2) = All_Candidates(2,i);
        M(2,3) = All_Candidates(1,i);
        M(2,4) = All_Candidates(6,i);
        M(3,5) = All_Candidates(3,i);
        M(3,6) = All_Candidates(4,i);
        M = M * dt;
        z = [x_truth(i+1);y_truth(i+1);z_truth(i+1)];
        s = [x_truth(i);y_truth(i);z_truth(i)];
        Sigma = Sigma + (z - M * Theta - s) * (z - M * Theta - s)';
    end
    Sigma = Sigma/(N-1).*eye(3);
    invSigma = inv(Sigma);
    D = zeros(6,6);
    c = zeros(6,1);
    for i = 1:N-1
            M = zeros(3,6);
        M(1,1) = All_Candidates(1,i);
        M(1,2) = All_Candidates(2,i);
        M(2,3) = All_Candidates(1,i);
        M(2,4) = All_Candidates(6,i);
        M(3,5) = All_Candidates(3,i);
        M(3,6) = All_Candidates(4,i);
        M = M * dt;
        z = [x_truth(i+1);y_truth(i+1);z_truth(i+1)];
        s = [x_truth(i);y_truth(i);z_truth(i)];
        D = D + M' * invSigma * M;
        c = c + M' * invSigma * (z - s);
    end
    D = D/(N-1);
    c = c/(N-1);
    Theta2 = D\c;
    lambda = (H/D*H')^(-1)*(H/D*c-g);
    Theta2_Constraint = D\(c-H'*lambda);
end
Sigma2 = Sigma; 

x2 = zeros(1,N);
y2 = zeros(1,N);
z2 = zeros(1,N);
x2(1) = 0;
y2(1) = 0;
z2(1) = 0;


rng(10)
for i = 2:N
    x2(i) = x2(i-1) + (Theta2_Constraint(1) * x2(i-1) + Theta2_Constraint(2) * y2(i-1)) * dt + sqrt(Sigma2(1,1)/dt) * sqrt(dt) * randn;
    y2(i) = y2(i-1) + (Theta2_Constraint(3) * x2(i-1) + Theta2_Constraint(4) * z2(i-1) * x2(i-1)) * dt + sqrt(Sigma2(2,2)/dt) * sqrt(dt) * randn;
    z2(i) = z2(i-1) + (Theta2_Constraint(5) * z2(i-1) + Theta2_Constraint(6) * x2(i-1) * y2(i-1)) * dt + sqrt(Sigma2(3,3)/dt) * sqrt(dt) * randn;
end
figure
subplot(3,5,1:3)
hold on
plot(dt:dt:N*dt, x2, 'g', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('x (identified model with a high threshold)','fontsize',12)
subplot(3,5,4)
ACF = autocorr(x2,Lag);
plot(0:dt:Lag*dt,ACF,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of x')
subplot(3,5,5)
[fi,xx] = ksdensity(x2);
plot(xx,fi,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of x')
subplot(3,5,[1:3]+5)
hold on
plot(dt:dt:N*dt, y2, 'g', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('y','fontsize',12)
subplot(3,5,4+5)
ACF = autocorr(y2,Lag);
plot(0:dt:Lag*dt,ACF,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of y')
subplot(3,5,5+5)
[fi,xx] = ksdensity(y2);
plot(xx,fi,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of y')
subplot(3,5,[1:3]+10)
hold on
plot(dt:dt:N*dt, z2, 'g', 'linewidth',2);
set(gca,'fontsize',12)
box on
title('z','fontsize',12)
xlabel('t')
subplot(3,5,4+10)
ACF = autocorr(z2,Lag);
plot(0:dt:Lag*dt,ACF,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('ACF of z')
subplot(3,5,5+10)
[fi,xx] = ksdensity(z2);
plot(xx,fi,'g','linewidth',2)
box on
set(gca,'fontsize',12)
title('PDF of z')