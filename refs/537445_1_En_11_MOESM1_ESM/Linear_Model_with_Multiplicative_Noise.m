% Chapter 7
% Linear model with multiplicative noise
% This code uses the linear model with multiplicative noise to characterize
% a one-dimensional observational time series, which is generated from a
% nonlinear cubic model

figure
% considering the two dynamical regimes
for j = 1:2
    if j == 1 % highly skewed regime
        a = -2; b = 2; c = 1; f = 0.1; A = 1; B = -1; sigma = 1; 
    elseif j == 2 % bimodal regime
        a = 2; b = 2; c = 1; f = 0.1; A = 1; B = -1; sigma = 1;
    end
    

    N = 300000;
    dt = 0.005;
    x = zeros(1,N);
    rng(11)% fix the random number seed
    for i = 2:N % generating the time series from a cubic nonlinear model
        x(i) = x(i-1) + (a * x(i-1) + b * x(i-1)^2 - c * x(i-1)^3 + f) * dt + ...
            (A - B * x(i-1)) * sqrt(dt) * randn + sigma * sqrt(dt) * randn;
    end
    [pdf_numerics, x_range] = ksdensity(x(1000:10:end)); % constructing the PDF

    mu = mean(x(1000:10:end)); % mean of the time series
    R = var(x(1000:10:end)); % variance of the time series
    pdf_G = normpdf(x_range, mu, sqrt(R)); % Gaussian fit of the PDF
    
    lag = 1500; % lag in computing the ACF
    ACF = autocorr(x,lag); % the ACF
    f = fit([0:dt:lag*dt]',ACF','exp1'); % using an exponential function to fit the ACF
    lambda = f.b; 
    Phi = zeros(1,length(x_range));
    for i = 2: length(x_range) % the function Phi in determining the multiplicative noise
        Phi(i) = trapz(x_range(1:i), (x_range(1:i) - mu) .* pdf_numerics(1:i));
    end
    sigma_mult = real( sqrt(2./pdf_numerics .* (lambda * Phi)) ); % the multiplicative noise coefficient
    rng(11)
    x2 = zeros(1,N);
    for i = 2:N % generate a time series from the calibrated linear model with multiplicative noise
        l = sum(x2(i-1) > x_range); % to deal with the boundary point issues
        if l < 1 
            l = 1;
        end
        if l > length(x_range)
            l = length(x_range);
        end
        sigma_x = sigma_mult(l);
        x2(i) = x2(i-1) + lambda * ( x2(i-1) - mu) * dt + sigma_x * sqrt(dt) * randn;
    end
    [pdf_numerics2, x_range] = ksdensity(x2(1000:10:end),x_range); % PDF of the simulated time series
    ACF2 = autocorr(x2,lag); % ACF of the simulated time series
    % compare the observed time series from the nonlinear model and the 
    % simulated time series from the linear model with multiplicative noise
    subplot(2,6,[1:3]+(j-1)*6)
    hold on
    plot(dt:dt:N*dt, x, 'b','linewidth',2)
    plot(dt:dt:N*dt, x2, 'r','linewidth',2)
    box on
    xlim([1000,1100])
    set(gca,'fontsize',16)
    if j == 1
        title('(a) Time series')
    end
    if j == 2
        xlabel('t')
    end
    % compare the ACFs and showing the exponential fit
    subplot(2,6,5+(j-1)*6)
    hold on
    plot(0:dt:lag*dt,ACF','b','linewidth',2);
    plot(0:dt:lag*dt,exp(lambda*[0:dt:lag*dt]),'k--','linewidth',2)
    plot(0:dt:lag*dt,ACF2','r','linewidth',2);
    box on
    set(gca,'fontsize',16)
    if j == 1
        title('(c) ACFs')
    end
    if j == 2
        xlabel('t (lag)')
    end
    % compare the PDFs and showing the Gaussian fit
    subplot(2,6,4+(j-1)*6)
    hold on
    plot(x_range, pdf_numerics, 'b' ,'linewidth',2)
    plot(x_range, pdf_G, 'k--','linewidth',2)
    plot(x_range, pdf_numerics2, 'r' ,'linewidth',2)
    box on
    set(gca,'fontsize',16)
    xlim([x_range(1),x_range(end)])
    if j == 2
        xlabel('x')
    end
    ylabel('p(x)')
    if j == 1
        title('(b) PDFs')
        legend('Cubic nonlinear model','Gaussian fit','Linear model with multiplicative noise')
    end
    % showing the computed multiplicative noise coefficient
    subplot(2,6,6+(j-1)*6)
    plot(x_range,sigma_mult,'r','linewidth',2)
    box on
    set(gca,'fontsize',16)
    if j == 1
        title('(d) \sigma(x)')
    end
    if j == 2
        xlabel('x')
    end
end