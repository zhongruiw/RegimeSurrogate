% Chapter 3
% Numerical integration of \int_0^1 x^2 dx
% The three methods described in the main text are provided here

% Method 1
NumRepeat = 1000; % repeat the experiments to quantify the uncertainty
disp('Method 1: see Figure 1')
figure(1)
variance_all = zeros(1,5); % save the variance of the estimator and study its relationship with N
for i = 1:5
    subplot(2,3,i)
    N = 10^i; % total number of samples 
    x = rand(NumRepeat,N);
    I = mean(x.^2,2); % the key step of the Monte Carlo simulation
    variance_all(i) = var(I);
    hold on
    % histogram is utilized to show the uncertainy with a fixed N exploiting the repeated experiments
    histogram(I,30,'facecolor','blue','edgecolor','blue') 
    plot([1/3,1/3],[0,100],'k--','linewidth',2)
    set(gca,'fontsize',12)
    xlim([0.1,0.7])
    box on
    title(['N = ', num2str(N)],'fontsize',16)
end
% This last subplot validates that the variance of the estimator is
% proportional to the inverse of N. In other words, the variance of the
% etimator times N should be a constant. 
N_all = 10.^[1:5];
subplot(2,3,6)
plot(N_all, variance_all .* N_all, '--ob', 'linewidth', 2)
box on
set(gca,'fontsize',12)
set(gca,'xscale','log')
xlabel('N')
ylabel('Variance of the estimator \times N')
title('Var \times N v.s. N')
ylim([0,0.5])
% Method 2
disp('Method 2:')
for i = 1:5
	N = 10^i; % total number of samples 
	x = rand(1,N); y = rand(1,N); % sample two-dimensional random numbers
	I = sum(x.^2 > y)/N; % count the percentange of the samples below the curve y = x^2
    disp(['N = ', num2str(N), '; I = ', num2str(I)]);
end

% Method 3
disp('Method 3: see Figure 3')
figure(3)
x = 0:0.01:1; y = x.^2;
plot(x,y,'b','linewidth',2)
grayColor = [.7 .7 .7];
hold on
NumBin = 10; % number of bin utilized 
for i = 1:NumBin
    xx = rand * 1/NumBin + 1/NumBin*(i-1); % sample a random number within each sub-interval
    yy = xx^2;
    patch([1/NumBin*(i-1), 1/NumBin*i, 1/NumBin*i, 1/NumBin*(i-1)], [0, 0, yy, yy], [.7 .7 .7], 'FaceAlpha', .5)
    plot([xx, xx], [0, yy], 'r', 'linewidth', 2)
end
