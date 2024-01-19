% Chapter 9
% Using MCMC (Metropolis algorithm) to sample a one-dimensional target 
% distribution

%%%%%%%% below is just plotting the figure with five points
x = -1.3:0.01:2.5; % range of the variable
f = exp(-3.2*x-11.2*x.^2+21.33*x.^3-8*x.^4); % function f that is proportional to the unknown probability p
p = f/trapz(x,f); % probability p
x_accept = [0,-0.1242,0.1530,-0.0432];
x_reject = [-1.1722,-0.6503];
p_accept = exp(-3.2*x_accept-11.2*x_accept.^2+21.33*x_accept.^3-8*x_accept.^4)/trapz(x,f);
p_reject = exp(-3.2*x_reject-11.2*x_reject.^2+21.33*x_reject.^3-8*x_reject.^4)/trapz(x,f);
figure
hold on
plot(x,p,'b','linewidth',2)
plot(x,f,'--k','linewidth',2)
plot(x_accept,p_accept,'ro','linewidth',2)
plot(x_reject,p_reject,'go','linewidth',2)
box on
set(gca,'fontsize',12)
xlim([x(1),x(end)])
legend('p(x)','f(x)','Points in the chain','Rejected points')
title('Metropolis algorithm in sampling a 1D distribution')
%%%%%%%%
% The actual MCMC start from here
rng(1) % fix the random number seed
K = 100000; % total number of samples
sample = zeros(1,K);
f0 = exp(-3.2*sample(1)-11.2*sample(1)^2+21.33*sample(1)^3-8*sample(1)^4); % evaluation at the initial point
for i = 2:K % MCMC
    sample(i) = sample(i-1) + randn;
    f1 = exp(-3.2*sample(i)-11.2*sample(i)^2+21.33*sample(i)^3-8*sample(i)^4); % evaluation
    rd = rand; % a uniform distributed random number
    if rd < f1/f0 % accept or reject
        f0 = f1;
    else
        sample(i) = sample(i-1);
    end
end
% plotting the PDF based on the constructed Markov chain
figure
for i = 1:4
    if i == 1
        K = 100;
    elseif i == 2
        K = 1000;
    elseif i == 3
        K = 10000;
    else
        K = 100000;
    end
    subplot(2,2,i)
    hold on
    plot(x,p,'b','linewidth',2)
    histogram(sample(1:K),30, 'Normalization','pdf');
    box on
    set(gca,'fontsize',16)
    xlim([x(1),x(end)])
    title(['N = ',num2str(K)]);
    if i == 1
        legend('Target','Reconstruction')
    end
end
 