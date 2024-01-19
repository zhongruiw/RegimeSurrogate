% Chapter 6
% Find the most sensitive direction via Fisher information and direct method 

% model parameters
% f = 1; a = 1; sigma = 1/2; % the first test
f = 1; a = 1; sigma = 5/2; % the second test
% mean and variance
mu = f/a;
C = sigma^2/2/a;

% Method I: using fisher information

MM = [1,-f/a;-f/a,f^2/a^2+sigma^2/4/a]; % the common factor a^3C has been removed
[eigvec,eigval] = eig(MM);
[num, index] = max(abs(diag(eigval))); % the most sensitive direction corresponding to the bigger eigenvalue 
disp('The most sensitive direction via Fisher information:')
disp(transpose(eigvec(:,index)));

% Method II: using direct numerical searching method
N = 601; % total number of grid points
delta = 0.05; % a small perturbation

% computing the difference using relative entropy
Error = zeros(N,N); % total difference measured by relative entropy
Error1 = zeros(N,N); % difference in the signal
Error2 = zeros(N,N); % difference in the dispersion
% information difference at all possible directions
for i = 1:N
    for j = 1:N
        dist = ((i - 301)/300)^2 + ((j-301)/300)^2;
        if dist<=1
            f_M = f+(i - 301)/300*delta;
            a_M = a+(j - 301)/300*delta;
            mu_M = f_M/a_M;
            C_M = sigma^2/2/a_M;
            Error(i,j) = 1/2 * (mu - mu_M)^2/C_M + 1/2 *(-log(C/C_M) + C/C_M - 1);
            Error1(i,j) = 1/2 * (mu - mu_M)^2/C_M ;
            Error2(i,j) = 1/2 *(-log(C/C_M) + C/C_M - 1);
        else
            Error(i,j) = nan;
            Error1(i,j) = nan;
            Error2(i,j) = nan;
        end
    end
end
% showing the results from Method II
figure
[xx,yy] = meshgrid([-1:1/300:1]*delta,[-1:1/300:1]*delta);
subplot(1,3,1) % total
contourf(xx,yy,Error',30,'linestyle','none');
[~,I] = max(Error(:));
[I_row, I_col] = ind2sub(size(Error),I);
Irow1 = (I_row-301)/300*delta; Irow2 = -Irow1;
Icol1 = (I_col-301)/300*delta; Icol2 = -Icol1;
hold on
plot([Irow1, Irow2], [Icol1, Icol2],'--k','linewidth',2)
disp('The most sensitive direction via direct searching:')
disp([Irow1,Icol1]/norm([Irow1,Icol1]));
h = colorbar;
cx = get(h,'Limits');
colormap jet
axis equal
title('Total','fontsize',18)
set(gca,'fontsize',12)
box on
xlabel('\delta{f}','fontsize',16)
ylabel('\delta{a}','fontsize',16)
subplot(1,3,2) % signal
contourf(xx,yy,Error1',30,'linestyle','none')
[~,I] = max(Error1(:));
[I_row, I_col] = ind2sub(size(Error),I);
Irow1 = (I_row-301)/300*delta; Irow2 = -Irow1;
Icol1 = (I_col-301)/300*delta; Icol2 = -Icol1;
hold on
plot([Irow1, Irow2], [Icol1, Icol2],'--k','linewidth',2)
colorbar
caxis([cx(1),cx(2)])
colormap jet
axis equal
title('Signal','fontsize',18)
set(gca,'fontsize',12)
box on
xlabel('\delta{f}','fontsize',16)
ylabel('\delta{a}','fontsize',16)
subplot(1,3,3) % dispersion
contourf(xx,yy,Error2',30,'linestyle','none')
[M,I] = max(Error2(:));
[I_row, I_col] = ind2sub(size(Error),I);
Irow1 = (I_row-301)/300*delta; Irow2 = -Irow1;
Icol1 = (I_col-301)/300*delta; Icol2 = -Icol1;
hold on
plot([Irow1, Irow2], [Icol1, Icol2],'--k','linewidth',2)
colorbar
caxis([cx(1),cx(2)])
colormap jet
axis equal
title('Dispersion','fontsize',18)
set(gca,'fontsize',12)
box on
xlabel('\delta{f}','fontsize',16)
ylabel('\delta{a}','fontsize',16)
