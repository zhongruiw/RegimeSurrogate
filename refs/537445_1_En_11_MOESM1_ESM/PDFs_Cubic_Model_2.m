% Chapter 4
% The equilibrium PDF of the cubic model
% dx = (f + ax + bx^2 - cx^3) dt + (A - Bx) dWC + sigma dWA
% Here A and B are zero
rng(1) % fix the random number seed
figure
for j = 1:4
    if j == 1 % nearly Gaussian regime
        a = -2.2; b = 0; c = 0; f = 2;  sigma = 1;
    elseif j == 2 % highly skewed regime
        a = -2; b = 2.5; c = 1; f = 0.0; sigma = 1;
    elseif j == 3 % sub-Gaussian regime
        a = -1; b = 0; c = 3; f = 0.0; sigma = 1;
    else % bimodal regime
        a = 2; b = -1; c = 2; f = 0.5; sigma = 1;
    end
     
    % A numerical simulation is also carried out here. The numerical
    % results can be shown along with the analytic solutions if needed.
    N = 300000; % total time points in numerical simulation 
    dt = 0.005; % numerical integration time step
    x = zeros(1,N); % state variable

    for i = 2:N % numerical integration
        x(i) = x(i-1) + (a * x(i-1) + b * x(i-1)^2 - c * x(i-1)^3 + f) * dt + ...
             + sigma * sqrt(dt) * randn;
    end
    [pdf_numerics, x_range] = ksdensity(x(1000:10:end)); % PDF from numerics
    % PDF from analytic solution
    pdf_theory = exp(2/sigma^2 * (f * x_range + a/2 * x_range.^2 + b/3 * x_range.^3 - c/4 * x_range.^4));
    pdf_theory = pdf_theory / trapz(x_range, pdf_theory);
    % computing the mean and variance to form the Gaussian fit
    mu = mean(x(1000:10:end));
    R = var(x(1000:10:end));
    pdf_G = normpdf(x_range, mu, sqrt(R)); % Gaussian fit of the PDF
    subplot(2,4,j)
    hold on
    plot(x_range, pdf_theory, 'b' ,'linewidth',2)
    plot(x_range, pdf_G, 'r--','linewidth',2)
    box on
    set(gca,'fontsize',16)
    xlim([x_range(1),x_range(end)])
    xlabel('x')
    ylabel('p(x)')
    if j == 1
        legend('Truth','Gaussian fit')
        title('(a) Nearly Gaussian Regime')
    elseif j == 2
        title('(b) Highly Skewed Regime')
    elseif j == 3
        title('Sub-Gaussian Regime')
    else
        title('Bimodal Regime')
    end
    subplot(2,4,4+j)% showing the PDF in the logarithm scale to better see the tail behavior 
    hold on
    plot(x_range, pdf_theory, 'b' ,'linewidth',2)
    plot(x_range, pdf_G, 'r--','linewidth',2)
    box on
    set(gca,'fontsize',16)
    set(gca,'yscale','log')
    ylim([1e-5,3])
    set(gca,'ytick',[1e-5,1])
    xlim([x_range(1),x_range(end)])
    xlabel('x')
    ylabel('p(x) in log scale')
end