% Chapter 6
% Quantifying predictability and model error using information criteria
% based on a linear Gaussian model
rng(11) % fix the random number seeds
t0 = 0; tf = 10; % initial time and final time
t = t0:0.1:tf; L = length(t); % time discretization for plotting the results
a1 = 0.5; f1 = 0; sigma1 = 1; % parameters in the linear Gaussian model
u0 = 3; % initial value of the truth
dt = 0.005; % numerical integration time step
N = round(tf/dt); % total number of points in time
u = zeros(1,N+1); % state variable
u(1) = u0; % initial value

for i = 2:N % numerical intergration of the linear Gaussian model
    u(i) = u(i-1) + (-a1 * u(i-1) + f1) * dt + sqrt(dt) * randn * sigma1;
end

figure
for k = 1:3 % computing the three cases with different sources of model error    
    if k == 1 % case 1
        a2 = 0.5; f2 = -0.5; sigma2 = 1;
    elseif k == 2 % case 2
        a2 = 1; f2 = 0; sigma2 = sqrt(2);
    else % case 3
        a2 = 0.5; f2 = 3; sigma2 = 1;
    end
    % the time evolution of the perfect system
    x0 = 3;
    x0v = 0.1;
    x1_mean = x0 * exp(-a1*t) + f1/a1 * (1-exp(-a1*t));
    x1_var = x0v * exp(-2*a1*t) + sigma1^2/2/a1 * (1-exp(-2*a1*t));
    x1_mean_eq = f1/a1;
    x1_var_eq = sigma1^2/2/a1;
    % the time evolution of the imperfect system
    x2_mean = x0 * exp(-a2*t) + f2/a2 * (1-exp(-a2*t));
    x2_var = x0v * exp(-2*a2*t) + sigma2^2/2/a2 * (1-exp(-2*a2*t));
    x2_mean_eq = f2/a2;
    x2_var_eq = sigma2^2/2/a2;

    Internal_1 = zeros(1,L); % internal predictability of the perfect system
    Internal_2 = zeros(1,L); % internal predictability of the imperfect system
    Model_Error = zeros(1,L); % model error between the two systems
    for i = 1:L
        Internal_1(i) = 1/2*(x1_mean(i)-x1_mean_eq)^2/x1_var_eq + 1/2 * (x1_var(i)/x1_var_eq - 1 - log(x1_var(i)/x1_var_eq));
        Internal_2(i) = 1/2*(x2_mean(i)-x2_mean_eq)^2/x2_var_eq + 1/2 * (x2_var(i)/x2_var_eq - 1 - log(x2_var(i)/x2_var_eq));
        Model_Error(i) = 1/2*(x1_mean(i)-x2_mean(i))^2/x2_var(i) + 1/2 * (x1_var(i)/x2_var(i) - 1 - log(x1_var(i)/x2_var(i)));
    end

    subplot(2,3,k)
    hold on
    graycolor = [0.7,0.7,0.7];
    x1_low = x1_mean-2*sqrt(x1_var);
    x1_high = x1_mean+2*sqrt(x1_var);
    x2_low = x2_mean-2*sqrt(x2_var);
    x2_high = x2_mean+2*sqrt(x2_var);
    h1 = patch([t,t(end:-1:1)],[x1_low,x1_high(end:-1:1)],graycolor,'edgecolor',graycolor,'facealpha',0.5,'edgealpha',0.5);
    h2 = patch([t,t(end:-1:1)],[x2_low,x2_high(end:-1:1)],'red','edgecolor','red','facealpha',0.3,'edgealpha',0.3);    
    h3 = plot(0:dt:N*dt, u,'color',graycolor,'linewidth',2);
    if k == 1        
        title('(a) a^M = 0.5, f^M = -0.5, \sigma^M = 1')
    elseif k == 2
        title('(b) a^M = 1, f^M = 0, \sigma^M = 2^{1/2}')
    else
        legend([h1,h2,h3],'Perfect model','Imperfect model','The signal true trajectory')
        title('(c) a^M = 0.5, f^M = 0, \sigma^M = 2')
    end
    ylim([-4.5,6])
    box on
    set(gca,'fontsize',16)
    subplot(2,3,k+3)
    yyaxis left
    hold on
    plot(t,Internal_1,'color',graycolor,'linewidth',2)
    h1 = plot(t,Internal_2,'-r','linewidth',2);
    ax = gca;
    ax.YColor = 'r';
    yyaxis right
    h2 = plot(t,Model_Error,'k','linewidth',2);
    ax = gca;
    ax.YColor = 'black';
    box on
    set(gca,'fontsize',16)
    xlabel('Lead time \tau')
    if k == 3
        legend('Internal predictability (perfect model)','Internal predictability (imperfect model)','Model error')
    end
end