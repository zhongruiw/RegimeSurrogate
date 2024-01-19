% Chapter 3
% Kernel density estimation of three different PDFs: Gaussian, fat-tailed,
% and bimodal

rng(1); % fix the random number generator to reproduce the results

% Test 1: Gaussian distribution
x = -5: 0.05: 5; % range of the domain
a = 0; b = 1; % parameters of mean and variance 
y = normpdf(x, a, b); % generate the corresponding Gaussian distribution
figure
for i = 1:4
    subplot(3,4,i) % different numbers of sample points
    if i == 1
        num = 2;
    elseif i == 2
        num = 4;
    elseif i == 3
        num = 10;
    else
        num = 1000;
    end
    z = randn(1,num);
    plot(x,y,'b','linewidth',2)
    [fi, xx, h] = ksdensity(z); % kernel density estimation
    hold on
    plot(xx, fi, 'r', 'linewidth',2)
    plot(z, ones(1,num) * 0.01, 'ko');
    if i <= 3
        for j = 1:num
            x0 = linspace(z(j) - 3*(h), z(j) + 3*(h), 100);
            y0 = 1 / num / h * 1/sqrt(2 * pi) * exp(- (x0 - z(j)).^2 / h^2); % for each KDE component
            plot(x0, y0, 'k', 'linewidth', 1)
        end
    end
    set(gca,'fontsize',16)
    box on
    title(['N = ', num2str(num)]);
    if i == 1
        legend('Truth', 'KDE approximation', 'Samples', 'Kernels')
    end
end

rng(1)% fix the random number generator to reproduce the results

% Test 2: a Gamma distribution, which has a fat tail
x = 0. :0.1: 25; % range of the domain
a = 2; b = 2; % parameters
y = gampdf(x, a, b); % generate the corresponding Gamma distribution

for i = 1:4
    subplot(3,4,4+i)% different numbers of sample points
    if i == 1
        num = 2;
    elseif i == 2
        num = 4;
    elseif i == 3
        num = 10;
    else
        num = 1000;
    end
    z = gamrnd(a, b, 1, num);
    plot(x, y, 'b', 'linewidth', 2)
    [fi, xx, h] = ksdensity(z); % kernel density estimation
    hold on
    plot(xx, fi, 'r', 'linewidth', 2)
    plot(z, ones(1,num) * 0.01, 'ko');
    if i <= 3
        for j = 1:num
            x0 = linspace(z(j) - 3*(h), z(j) + 3*(h), 100);
            y0 = 1 / num / h * 1/sqrt(2 * pi) * exp(- (x0 - z(j)).^2 / h^2); % for each KDE component
            plot(x0, y0, 'k', 'linewidth', 1)
        end
    end
    set(gca,'fontsize',16)
    box on
end
 

rng(1) % fix the random number generator to reproduce the results

% Test 3: bimodal distribution
x = -6: 0.05: 6; % range of the domain
gm = gmdistribution([-2, 2], [0.5, 0.5]); % build the bimodal distribution via a two-component Gaussian mixture
temp = random(gm, 10000); temp = reshape(temp, 1, []); [y,temp0] = ksdensity(temp, x);
num = 200;
z = random(gm, num/2); z = reshape(z, 1, []); % generate the truth, which is based on a large number (10000) of sample points
subplot(3, 4, 8 + 2)
plot(x,y,'b','linewidth',2)
[fi,xx, h] = ksdensity(z); % standard kernel denstify estimation
hold on
plot(xx,fi,'r','linewidth',2)
[bandwidth, density, xmesh, cdf] = kde(z, num, -6, 6); % solve-the-equation plug-in method
plot(xmesh, density, 'g', 'linewidth',2)
box on
set(gca,'fontsize',16)
legend('Truth', 'KDE with rule of thumb', 'KDE with solve-the-equation plug-in')

subplot(3, 4, 8 + 3); % compare the kernel components in the two methods
hold on
x0 = linspace(0 - 3*(h), 0 + 3*(h), 100);
y0 =  1/sqrt(2 * pi) * exp(-(x0 - 0).^2 / h^2);
plot(x0,y0,'r','linewidth',1)
x0 = linspace(0 - 3 * bandwidth, 0 + 3 * bandwidth, 100);
y0 = 1/sqrt(2 * pi) * exp(-(x0 - 0).^2 / bandwidth^2);
plot(x0, y0, 'g', 'linewidth', 1)
box on
set(gca,'fontsize',16);