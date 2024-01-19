% Chapter 8
% Using the hybrid algorithm to solve both the transient and equilibrium
% PDFs (namely solving the Fokker-Planck equation) 
% The test model here is the noisy Lorenz 63 model

rng(11); % fix the random number seed
% parameters in the noisy Lorenz 63 model
sg = 10;  rho = 28; beta = 8/3; % parameters in the deterministic part
sgm1 = 20; sgm2 = 20; sgm3 = 20; % parameters in the noise part

dt = 0.005; % numerical integration time step
N = 70; % total number of steps (transient phase)
% N = 1000; % total number of steps (equilibrium phase)
time_instant = N;

Num_MC = 150000; % number of samples in the Monte Carlo simulation
% initial values
x1_MC0 = randn(Num_MC,1)/2; 
x2_MC0 = randn(Num_MC,1)/2; 
x3_MC0 = randn(Num_MC,1)/2; 

% number of samples in the hybrid algorithm
Num = 500;
% MM is the numerical of points which is used only to plot the figure
MM = 150000;

% generating the true solution using Monte Carlo simulation
j = 2;
for i = 2: N
    x1_MC1 = x1_MC0 + sg * ( x2_MC0 - x1_MC0 ) * dt + sgm1 * randn(Num_MC,1) * sqrt(dt);
    x2_MC1 = x2_MC0 + ( x1_MC0 .* ( rho - x3_MC0 ) - x2_MC0 ) * dt + sgm2 * randn(Num_MC,1) * sqrt(dt);
    x3_MC1 = x3_MC0 + ( x1_MC0 .* x2_MC0 - beta * x3_MC0) * dt + sgm3 * randn(Num_MC,1) * sqrt(dt);
    if i == time_instant
        x1_save = x1_MC1;
        x2_save = x2_MC1;
        x3_save = x3_MC1;
    end
    x1_MC0 = x1_MC1;  x2_MC0 = x2_MC1;  x3_MC0 = x3_MC1; 
    j = j + 1;
end

% find the range of each variable 
[fi1t,xx1] = ksdensity(x1_save);
[fi2t,xx2] = ksdensity(x2_save);
[fi3t,xx3] = ksdensity(x3_save);
[fi1t,xx1] = ksdensity(x1_save,linspace(xx1(1),xx1(end),2^7));

[bandwidth1,fi1t,xmesh1,cdf1]=kde(x1_save,2^7,min(xx1),max(xx1)); % KDE using solve-the-equation plug-in method for x



figure
% showing two-dimensional PDFs of the truth (i.e., from Monte Carlo
% simulation)
subplot(2,6,10)
x = x1_save;
y = x2_save;
x(x<xx1(1)) = xx1(1);
x(x>xx1(end)) = xx1(end);
y(y<xx2(1)) = xx2(1);
y(y>xx2(end)) = xx2(end);
xi = xx1;
yi = xx2; 
xr = interp1(xi, 0.5:numel(xi)-0.5, x, 'nearest');        
yr = interp1(yi, 0.5:numel(yi)-0.5, y, 'nearest');
Z = accumarray([yr xr] + 0.5, 1, [length(yi) length(xi)]);
zz = Z/MM/(xx1(2)-xx1(1))/(xx2(2)-xx2(1));
imagesc(xi,yi,zz)
set(gca,'fontsize',12)
xlabel('x','FontName','Times','fontsize',16);
ylabel('y','FontName','Times','fontsize',16);
title('Truth p(x,y)','FontName','Times','fontsize',16);
zz(zz<1e-5) = 1e-5; fi12t = zz/trapz(xx1,trapz(xx2,zz));

subplot(2,6,11)
x = x2_save;
y = x3_save;
x(x<xx2(1)) = xx2(1);
x(x>xx2(end)) = xx2(end);
y(y<xx3(1)) = xx3(1);
y(y>xx3(end)) = xx3(end);
xi = xx2;
yi = xx3; 
xr = interp1(xi, 0.5:numel(xi)-0.5, x, 'nearest');        
yr = interp1(yi, 0.5:numel(yi)-0.5, y, 'nearest');
Z = accumarray([yr xr] + 0.5, 1, [length(yi) length(xi)]);
zz = Z/MM/(xx2(2)-xx2(1))/(xx3(2)-xx3(1));
imagesc(xi,yi,zz)
set(gca,'fontsize',12)
xlabel('y','FontName','Times','fontsize',16);
ylabel('z','FontName','Times','fontsize',16);
title('Truth p(y,z)','FontName','Times','fontsize',16);
zz(zz<1e-5) = 1e-5; fi23t = zz/trapz(xx2,trapz(xx3,zz));

subplot(2,6,12)
x = x3_save;
y = x1_save;
x(x<xx3(1)) = xx3(1);
x(x>xx3(end)) = xx3(end);
y(y<xx1(1)) = xx1(1);
y(y>xx1(end)) = xx1(end);
xi = xx3;
yi = xx1; 
xr = interp1(xi, 0.5:numel(xi)-0.5, x, 'nearest');        
yr = interp1(yi, 0.5:numel(yi)-0.5, y, 'nearest');
Z = accumarray([yr xr] + 0.5, 1, [length(yi) length(xi)]);
zz = Z/MM/(xx3(2)-xx3(1))/(xx1(2)-xx1(1));
imagesc(xi,yi,zz)
set(gca,'fontsize',12)
xlabel('z','FontName','Times','fontsize',16);
ylabel('x','FontName','Times','fontsize',16);
title('Truth p(z,x)','FontName','Times','fontsize',16);
zz(zz<1e-5) = 1e-5; fi31t = zz/trapz(xx3,trapz(xx1,zz));


% nn is the number of points used to create a mixture component when
% showing the figure. Again, nn has nothing to do with the algorithm. It is
% only used for plotting the figure.
nn = round(MM/Num);

% below is the hybrid algorithm

x1_truth = zeros(Num, N);
x2_truth = zeros(Num, N);
x3_truth = zeros(Num, N);

% initial value
x1_truth0 = randn(Num,1)/2;
x2_truth0 = randn(Num,1)/2;
x3_truth0 = randn(Num,1)/2;

x1_truth(:,1) = x1_truth0;
x2_truth(:,1) = x2_truth0;
x3_truth(:,1) = x3_truth0;


% running the model numerically forward but only with a small number of
% ensemble Num.
j = 2;
for i = 2: N

    x1_truth1 = x1_truth0 + sg * ( x2_truth0 - x1_truth0 ) * dt + sgm1 * randn(Num,1) * sqrt(dt);
    x2_truth1 = x2_truth0 + ( x1_truth0 .* ( rho - x3_truth0 ) - x2_truth0 ) * dt + sgm2 * randn(Num,1) * sqrt(dt);
    x3_truth1 = x3_truth0 + ( x1_truth0 .* x2_truth0 - beta * x3_truth0 ) * dt + sgm3 * randn(Num,1) * sqrt(dt);
    

    x1_truth0 = x1_truth1;  x2_truth0 = x2_truth1;  x3_truth0 = x3_truth1; 




    x1_truth(:,j) = x1_truth1;
    x2_truth(:,j) = x2_truth1;
    x3_truth(:,j) = x3_truth1;
    j = j + 1;

end

x1_trace = x1_truth;
% using KDE for the observed variable x
[bandwidth,fi1,xmesh1,cdf] = kde(x1_truth(:,time_instant),2^7,min(xx1),max(xx1));
 

gamma_mean_trace = zeros(2*Num,N); % Store the posterior mean at each time step
gamma_cov_trace = zeros(4*Num,N); % Store the posterior cov at each time step
gamma_mean0 = zeros(2*Num, 1); % posterior mean at previous step
gamma_cov0 = 2*eye(2*Num); % posterior cov at previous step; only the block diagonal 2*2 pieces are useful
mark1 = 1: 2: 2*Num-1; % Index just for the convenience of saving each posterior matrix 
mark2 = 2: 2: 2*Num; % Same as mark1
gamma_cov_trace(1:4:end,1) = diag(gamma_cov0(mark1,mark1)); % Store the posterior matrices; 
gamma_cov_trace(2:4:end,1) = diag(gamma_cov0(mark1,mark2));
gamma_cov_trace(3:4:end,1) = diag(gamma_cov0(mark2,mark1));
gamma_cov_trace(4:4:end,1) = diag(gamma_cov0(mark2,mark2));

% noise matrix; b1 is Sigma_Y; InvBoB is the inverse of Sigma_X * Sigma_X^-1; 
% putting all Num ensembles into a big matrix to avoid the for loop
b1 = zeros(2*Num,1);
b1(1:2:end) = sgm2;
b1(2:2:end) = sgm3;
b1 = diag(b1);
invBoB = 1 / sgm1^2 * eye(Num);

%%% Constant matrix a0 
a0 = zeros(2*Num,1);
% x0 is x at previous step
x0 = x1_trace(:,1);
% solving the conditional distribution using the filtering method
for s = 2:N
    if mod(s,400)==1
        disp(s)
    end
    % x is u1 at current step; x and x0 are needed to update posterior mean
    x = x1_trace(:,s);
    
    a1 = zeros(2*Num,1);
    a1(1:2:end) = -1;
    a1(2:2:end) = -beta;
    a1 = diag(a1);
    for i = 1:2:2*Num-1
        a1(i,i+1) = - x1_trace((i+1)/2,s-1);
        a1(i+1,i) = x1_trace((i+1)/2,s-1);
    end
    % A0 and A1 are in uI equation; special trick to put each of both into 
    % a big matrix
    A0 = -sg * x1_trace(:,s-1);
    A1 = zeros(2*Num^2,1);
    A1(1:2*(Num+1):end) = sg;
    A1 = transpose(reshape(A1,2*Num,Num));
    % update a0 at each step, since a0 depends on uI
    a0(1:2:end) = rho * x1_trace(:,s-1);
    
    % update posterior mean and cov; eqn 3 in CM paper
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (x-x0 - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
    
    % Store the poseterior mean and cov at each step
    gamma_mean_trace(:,s) = gamma_mean;
    gamma_cov_trace(1:4:end,s) = diag(gamma_cov(mark1,mark1));
    gamma_cov_trace(2:4:end,s) = diag(gamma_cov(mark1,mark2));
    gamma_cov_trace(3:4:end,s) = diag(gamma_cov(mark2,mark1));
    gamma_cov_trace(4:4:end,s) = diag(gamma_cov(mark2,mark2));
    
    % Update for the next iteration step
    x0 = x;
    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;
end


% plotting the results from the hybrid method
% plotting 2D distributions
subplot(2,6,5)
% the conditional Gaussian mixture
MU = reshape(gamma_mean_trace(:,time_instant),2,[])'; % posterior mean
SIGMA = zeros(2,2,Num);
for i = 1:Num % posterior covariance
    SIGMA(:,:,i) = reshape(gamma_cov_trace(4*i-3:4*i,time_instant),2,2)';
end
p = ones(1,Num)/Num; % weight
% Conditional Gaussian mixture
obj = gmdistribution(MU,SIGMA,p);
Y0 = random(obj,MM);
x = Y0(:,1);
y = Y0(:,2);
x(x<xx2(1)) = xx2(1);
x(x>xx2(end)) = xx2(end);
y(y<xx3(1)) = xx3(1);
y(y>xx3(end)) = xx3(end);
xi = xx2;
yi = xx3; 
xr = interp1(xi, 0.5:numel(xi)-0.5, x, 'nearest');        
yr = interp1(yi, 0.5:numel(yi)-0.5, y, 'nearest');
Z = accumarray([yr xr] + 0.5, 1, [length(yi) length(xi)]);
zz = Z/MM/(xx2(2)-xx2(1))/(xx3(2)-xx3(1));
imagesc(xi,yi,zz)
set(gca,'fontsize',12)
xlabel('y','FontName','Times','fontsize',16);
ylabel('z','FontName','Times','fontsize',16);
title('Recovered p(y,z)','FontName','Times','fontsize',16);
zz(zz<1e-5) = 1e-5; fi23 = zz/trapz(xx2,trapz(xx3,zz));
[fi2,xx2] = ksdensity(Y0(:,1),xx2);
[fi3,xx3] = ksdensity(Y0(:,2),xx3);

subplot(2,6,4)
MU = reshape(gamma_mean_trace(:,time_instant),2,[])';
MU_temp = MU;
MU_temp(:,2) = MU(:,1);
MU_temp(:,1) = x1_truth(:,time_instant)';
MU = MU_temp;
for j = 1:Num
    SIGMA = reshape(gamma_cov_trace(4*j-3:4*j,time_instant),2,2)';
    SIGMA_temp = SIGMA*0;
    SIGMA_temp(2,2) = SIGMA(1,1);
    SIGMA_temp(1,1) = bandwidth(1)^2;
    SIGMA = SIGMA_temp;
    Y0((j-1)*nn+1:j*nn,:) = mvnrnd(MU(j,:),SIGMA,nn);
end
x = Y0(:,1);
y = Y0(:,2);
x(x<xx1(1)) = xx1(1);
x(x>xx1(end)) = xx1(end);
y(y<xx2(1)) = xx2(1);
y(y>xx2(end)) = xx2(end);
xi = xx1;
yi = xx2; 
xr = interp1(xi, 0.5:numel(xi)-0.5, x, 'nearest');        
yr = interp1(yi, 0.5:numel(yi)-0.5, y, 'nearest');
Z = accumarray([yr xr] + 0.5, 1, [length(yi) length(xi)]);
zz = Z/MM/(xx1(2)-xx1(1))/(xx2(2)-xx2(1));
imagesc(xi,yi,zz)
set(gca,'fontsize',12)
xlabel('x','FontName','Times','fontsize',16);
ylabel('y','FontName','Times','fontsize',16)
title('Recovered p(x,y)','FontName','Times','fontsize',16);
zz(zz<1e-5) = 1e-5; fi12 = zz/trapz(xx1,trapz(xx2,zz));

subplot(2,6,6)
MU = reshape(gamma_mean_trace(:,time_instant),2,[])';
MU_temp = MU;
MU_temp(:,2) = MU(:,2);
MU_temp(:,1) = x1_truth(:,time_instant)';
MU = MU_temp;
for j = 1:Num
    SIGMA = reshape(gamma_cov_trace(4*j-3:4*j,time_instant),2,2)';
    SIGMA_temp = SIGMA*0;
    SIGMA_temp(2,2) = SIGMA(2,2);
    SIGMA_temp(1,1) = bandwidth(1)^2;
    SIGMA = SIGMA_temp;
    Y0((j-1)*nn+1:j*nn,:) = mvnrnd(MU(j,:),SIGMA,nn);

end
x = Y0(:,2);
y = Y0(:,1);
x(x<xx3(1)) = xx3(1);
x(x>xx3(end)) = xx3(end);
y(y<xx1(1)) = xx1(1);
y(y>xx1(end)) = xx1(end);
xi = xx3;
yi = xx1; 
xr = interp1(xi, 0.5:numel(xi)-0.5, x, 'nearest');        
yr = interp1(yi, 0.5:numel(yi)-0.5, y, 'nearest');
Z = accumarray([yr xr] + 0.5, 1, [length(yi) length(xi)]);
zz = Z/MM/(xx1(2)-xx1(1))/(xx3(2)-xx3(1));
imagesc(xi,yi,zz)
set(gca,'fontsize',12)
xlabel('z','FontName','Times','fontsize',16);
ylabel('x','FontName','Times','fontsize',16)
title('Recovered p(z,x)','FontName','Times','fontsize',16);
zz(zz<1e-5) = 1e-5; fi31 = zz/trapz(xx3,trapz(xx1,zz));

% plotting 1D distribution
subplot(2,6,1)
hold on
plot(xx1,fi1,'b','linewidth',2)
plot(xx1,fi1t,'r--','linewidth',2)
set(gca,'fontsize',12)
box on
title('x','FontName','Times','fontsize',16);
legend('Recovered','Truth')
subplot(2,6,2)
hold on
plot(xx2,fi2,'b','linewidth',2)
plot(xx2,fi2t,'r--','linewidth',2)
set(gca,'fontsize',12)
box on
title('y','FontName','Times','fontsize',16);
subplot(2,6,3)
hold on
plot(xx3,fi3,'b','linewidth',2)
plot(xx3,fi3t,'r--','linewidth',2)
set(gca,'fontsize',12)
box on
title('z','FontName','Times','fontsize',16);

% computing the error in the 1D and 2D distributions using relative entropy 
modelerror1 = trapz(xx1, fi1t.*log(fi1t./fi1));
modelerror2 = trapz(xx2, fi2t.*log(fi2t./fi2));
modelerror3 = trapz(xx3, fi3t.*log(fi3t./fi3));
modelerror4 = trapz(xx1, trapz(xx2,fi12t.*log(fi12t./fi12)));
modelerror5 = trapz(xx2, trapz(xx3,fi23t.*log(fi23t./fi23)));
modelerror6 = trapz(xx3, trapz(xx1,fi31t.*log(fi31t./fi31)));
disp('Error (via relative entropy) in   p(x),   p(y),   p(z),   p(x,y),   p(y,z),   p(z,x)')
disp([modelerror1, modelerror2, modelerror3, modelerror4, modelerror5, modelerror6])

 