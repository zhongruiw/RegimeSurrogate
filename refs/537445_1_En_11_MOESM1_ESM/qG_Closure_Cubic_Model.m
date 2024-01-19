% Chapter 4
% quasi-Gaussian closure for the moment equations of the cubic model
% dx = (f + ax + bx^2 - cx^3) dt + (A - Bx) dWC + sigma dWA

rng(1) % fix the random number seed
figure
Ens = 1500; % ensemble size in the Monte Carlo simulation for comparison
N = 4000; % total number of time steps 
dt = 0.005; % numerical integration time step

x_range = linspace(-3,6,200); % the range of the state variable for showing the PDF
gap = 10; % record the solution at every 'gap' steps 
xpdf_all = zeros(200,round(N/gap)); % time evolution of the PDF


% showing the results in different dynamical regimes
for j = 1:4 
    if j == 1 % nearly Gaussian regime
        a = -2.2; b = 0; c = 0; f = 2; A = 0.1; B = 0.1; sigma = 1;
    elseif j == 2 % highly skewed regime
        a = -4; b = 2; c = 1; f = 0.1; A = 1; B = -1; sigma = 1;
    elseif j == 3 % fat-tailed regime
        a = -3; b = -1.5; c = 0.5; f = 0.0; A = 0.5; B = -1; sigma = 1;
    else % bimodal regime
        a = 4; b = 2; c = 1; f = 0.1; A = 1; B = -1; sigma = 1;
    end
    
    % solution from the Monte Carlo simulation
    x_all = zeros(Ens,round(N/gap));
    x_old = randn(Ens,1) * 0 - 2;
    x_all(:,1) = x_old;
    % qG closure with the same initial value
    x_mean_all = zeros(1,round(N/gap));
    x_var_all = zeros(1,round(N/gap));    
    x_mean_old = 0 - 2;
    x_var_old = 0;
    x_mean_all(1) = x_mean_old;
    x_var_all(1) = x_var_old;
    % qG closure with the a different initial value (to check the long-term
    % behavior)
    x2_mean_all = zeros(1,round(N/gap));
    x2_var_all = zeros(1,round(N/gap));
    x2_mean_old = 0 + 2;
    x2_var_old = 1;
    x2_mean_all(1) = x2_mean_old;
    x2_var_all(1) = x2_var_old;
    % bare truncation solution of the moment equations (with the same
    % initial condition as the truth)
    x3_mean_all = zeros(1,round(N/gap));
    x3_var_all = zeros(1,round(N/gap));
    x3_mean_old = 0 - 2;
    x3_var_old = 0;
    x3_mean_all(1) = x3_mean_old;
    x3_var_all(1) = x3_var_old;
    k = 2;
    for i = 2:N % time evolution
        % Monte Carlo simulation
        x_new = x_old + (a * x_old + b * x_old.^2 - c * x_old.^3 + f) * dt + ...
            (A - B * x_old) * sqrt(dt) .* randn(Ens,1) + sigma * sqrt(dt) * randn(Ens,1);
        % qG closure with the same initial value
        x_mean_new = x_mean_old + (a * x_mean_old + b * x_mean_old^2 + b * x_var_old - c * x_mean_old^3 ...
            - 3 * c * x_mean_old * x_var_old + f) * dt; 
        x_var_new = x_var_old + ( (2 * ( a + 2 * b * x_mean_old - 3 * c * x_mean_old^2 - 3 * c * x_var_old) + B^2) * x_var_old ...
            + ( A^2 + B^2 * x_mean_old^2 - 2 * A * B * x_mean_old + sigma^2) ) * dt;
        x_old = x_new;
        x_mean_old = x_mean_new;
        x_var_old = x_var_new;
        % qG closure with the a different initial value
        x2_mean_new = x2_mean_old + (a * x2_mean_old + b * x2_mean_old^2 + b * x2_var_old - c * x2_mean_old^3 ...
            - 3 * c * x2_mean_old * x2_var_old + f) * dt; 
        x2_var_new = x2_var_old + ( (2 * ( a + 2 * b * x2_mean_old - 3 * c * x2_mean_old^2 - 3 * c * x2_var_old) + B^2) * x2_var_old ...
            + ( A^2 + B^2 * x2_mean_old^2 - 2 * A * B * x2_mean_old + sigma^2) ) * dt;
        x2_mean_old = x2_mean_new;
        x2_var_old = x2_var_new;
        % bare truncation solution of the moment equations
        x3_mean_new = x3_mean_old + (a * x3_mean_old + b * x3_mean_old^2 + b * x3_var_old - c * x3_mean_old^3 ...
            - 3 * c * x3_mean_old * x3_var_old + f) * dt; 
        x3_var_new = x3_var_old + ( (2 * ( a + 2 * b * x3_mean_old - 3 * c * x3_mean_old^2 ) + B^2) * x3_var_old ...
            + ( A^2 + B^2 * x3_mean_old^2 - 2 * A * B * x3_mean_old + sigma^2) ) * dt;
        x3_mean_old = x3_mean_new;
        x3_var_old = x3_var_new;
        % Record the solution at every 'gap' steps 
        if mod(i,gap) == 1
            x_all(:,k) = x_old;  
            x_mean_all(k) = x_mean_old;
            x_var_all(k) = x_var_old;
            x2_mean_all(k) = x2_mean_old;
            x2_var_all(k) = x2_var_old;
            x3_mean_all(k) = x3_mean_old;
            x3_var_all(k) = x3_var_old;
            [fi,x_range] = ksdensity(x_old,x_range);
            xpdf_all(:,k) = fi';
            k = k + 1;
        end
    end

    grayColor = [.7 .7 .7];
    % the time evolution of the PDF
    subplot(3,4,j)
    [x1,y1] = meshgrid(x_range,dt:gap*dt:N*dt);
    mesh(x1,y1,xpdf_all')
    xlim([-3,6])
    set(gca,'fontsize',16)
    box on
    if j == 1
        title('(a) Nearly Gaussian Regime')
    elseif j == 2
        title('(b) Highly Skewed Regime')
    elseif j == 3
        title('(c) Fat-tailed Regime')
    else
        title('(d) Bimodal Regime')
    end
    xlabel('x')
    ylabel('t')
    zlabel('p(x,t)');
    % time evolution of the moments using different methods
    subplot(3,4,j+4)
    hold on
    plot(dt:gap*dt:N*dt, mean(x_all),'b','linewidth',2)
    plot(dt:gap*dt:N*dt, x_mean_all,'r','linewidth',2)
    plot(dt:gap*dt:N*dt, x2_mean_all,'g','linewidth',2)
    plot(dt:gap*dt:N*dt, x3_mean_all,'--k','linewidth',2)
    if j <= 3
        xlim([0,4]);
    end
    if j == 1
        legend('Truth via Monte Carlo','qG closure with the true x(0)','qG closure with a differnt x(0)','Bare truncation')
    end
    box on
    set(gca,'fontsize',16)
    % time evolution of the variance using different methods
    subplot(3,4,j+8)
    hold on
    plot(dt:gap*dt:N*dt, var(x_all),'b','linewidth',2)
    plot(dt:gap*dt:N*dt, x_var_all,'r','linewidth',2)
    plot(dt:gap*dt:N*dt, x2_var_all,'g','linewidth',2)
    plot(dt:gap*dt:N*dt, x3_var_all,'--k','linewidth',2)
    if j <= 3
        xlim([0,4]);
    end
    box on
    set(gca,'fontsize',16)
    xlabel('t')
end