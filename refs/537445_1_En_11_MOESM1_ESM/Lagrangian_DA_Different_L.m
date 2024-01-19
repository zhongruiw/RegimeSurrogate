% Chapter 8
% Lagrangian data assimilation for testing the error and uncertainty
% decaying as a function f the number of tracers L
% This code should be run together with Lagrangian_DA.m, after running the
% Shallow_Water_Equation.m
figure
RE_s = zeros(1,7); % signal part of relative entropy 
RE_d = zeros(1,7); % dispersion part of relative entropy
LL = [2,5,10,20,30,50,100];  % number of tracers L in the test
for jj = 1:7
    disp(jj);
    L = LL(jj);
    Lagrangian_DA % call the code for running Lagrangian data assimilation with a given L
    RE_s(jj) = Relative_Entropy_Signal_All; % relative entropy in signal part
    RE_d(jj) = Relative_Entropy_Dispersion_All; % relative entropy in dispersion part
    if jj == 1 || 3 || 5 || 7 % only showing the results with every other test
        subplot(2,4,(jj+1)/2) % showing GB modes
        hold on
        plot(dt:dt:N*dt, real(u_hat(Dim_Ug*2+6,:)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean(Dim_Ug*2+6,:)), 'r', 'linewidth',2)
        patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(Dim_Ug*2+6,:))+2*sqrt(real(u_post_cov(Dim_Ug*2+6,:))), real(u_post_mean(Dim_Ug*2+6,end:-1:1))-2*sqrt(real(u_post_cov(Dim_Ug*2+6,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
        title(['L = ', num2str(L), ';  GB mode ( ', num2str(kk(1,6)),' , ', num2str(kk(2,6)), ' )'],'fontsize',16)
        set(gca,'fontsize',12)
        box on

        subplot(2,4,4+(jj+1)/2) % showing gravity modes
        hold on
        plot(dt:dt:N*dt, real(u_hat(6,:)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean(6,:)), 'r', 'linewidth',2)
        patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(6,:))+2*sqrt(real(u_post_cov(6,:))), real(u_post_mean(6,end:-1:1))-2*sqrt(real(u_post_cov(6,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
        title(['Gravity mode ( ', num2str(kk(1,6)),' , ', num2str(kk(2,6)), ' )'],'fontsize',14)
        set(gca,'fontsize',12)
        box on
        xlabel('t')
    end
    
end
% showing the information gain (uncertainty reduction) as a function of L
figure
hold on
plot(LL,RE_s,'b','linewidth',2)
plot(LL,RE_d,'r','linewidth',2)
legend('Signal','Dispersion')
set(gca,'fontsize',16)
box on
xlabel('L')
ylabel('Bits')
title('Information gain')