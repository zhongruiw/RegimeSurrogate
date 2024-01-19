% Chapter 1
% Simulation two-state Markov jump process
% The two states are named as 'st' and 'un'

rng(1); % fix random number seeds to reproduce results 
N = 250000; % total length of the simulation points
dt = 0.005; % numerical integration time step
X = zeros(1,N); % the simulated time series

mu = 0.5; % transition probability from 'un' to 'st'
nu = 0.1; % transition probability from 'st' to 'un'
X_st = 2.27; % state value at the 'st' state
X_un = -0.04; % value at the 'un' state

X(1) = X_st; % initial value

% simulate the time series
for i = 2:N
    rd = rand; % generate a random number to determine if the transition occurs
    if X(i-1) == X_un % current: 'un' state
        if rd < mu * dt 
            X(i) = X_st; % transition to 'st' state
        else
            X(i) = X_un; % remain as 'un'
        end
    else % current: 'st' state
        if rd < nu * dt 
            X(i) = X_un; % transition to 'un' state
        else
            X(i) = X_st; % remain as 'st'
        end
    end
end
% plot the results
figure
subplot(2,1,1)
plot(dt:dt:N*dt, X, 'b', 'linewidth',2)
ylim([X_un-0.1,X_st+0.1])
xlim([0,N*dt])
xlabel('t')
ylabel('States')
box on
set(gca,'fontsize',12)
title('Simulated time series')

% Compute the statistics

% Equilibrium distribution
p_eq_un_Numerical = length(find(X == X_un))/N; % equilibrium probability at 'un' state
p_eq_st_Numerical = length(find(X == X_st))/N; % equilibrium probability at 'st' state
p_eq_un_Theoretical = nu / (nu + mu); % equilibrium probability at 'un' state
p_eq_st_Theoretical = mu / (nu + mu); % equilibrium probability at 'st' state
disp('Equilibrium distribution:')
disp(['p_{eq}(un) theoretical result: ', num2str(p_eq_un_Theoretical), '; p_{eq}(st) theoretical result: ', num2str(p_eq_st_Theoretical), '.'])
disp(['p_{eq}(un) numerical result: ', num2str(p_eq_un_Numerical), '; p_{eq}(st) numerical result: ', num2str(p_eq_st_Numerical), '.'])

% Expectation
Mean_X_Numerical = mean(X); % taking the mean from the simulation
Mean_X_Theoretical = (mu * X_st + nu * X_un) / (mu + nu); % theoretical result
disp('Expectation:')
disp(['E[X] numerical result: ', num2str(Mean_X_Numerical), '; E[X] theoretical result: ', num2str(Mean_X_Theoretical), '.'])

% Switching times
T_st_Numerical = []; % save the switching time
T_un_Numerical = []; % save the switching time

% switching to the un state
flag = 0;
for i = 1:N
    if X(i) == X_st && flag == 0 % current state
        j = i; flag = 1;
    end
    if X(i) == X_un && flag == 1 % switching occurs
        T_st_Numerical = [T_st_Numerical, i-j]; % add this record to the switching time 
        flag = 0;
    end
end

% switching to the st state
flag = 0;
for i = 1:N
    if X(i) == X_un && flag == 0 % current state
        j = i; flag = 1;
    end
    if X(i) == X_st && flag == 1 % switching occurs
        T_un_Numerical = [T_un_Numerical, i-j]; % add this record to the switching time 
        flag = 0;
    end
end

T_st_Numerical = T_st_Numerical * dt; % multiplying by the time step
T_un_Numerical = T_un_Numerical * dt; % multiplying by the time step

disp('Switching times:')
disp('See the figure')

subplot(2,2,3)
hold on
histogram(T_st_Numerical, 20, 'EdgeColor','b','FaceColor','b','normalization', 'pdf'); % numerical result
t_range = dt: dt: max(T_st_Numerical); 
plot(t_range, exp(- nu * t_range) * nu, 'r', 'linewidth', 2); % theoretical result; multiplying nu to make it to be a PDF
box on
set(gca,'fontsize',12)
title('Switching time from st to un states')
xlabel('t')
ylabel('Probability')
legend('Numerical','Theoretical')

subplot(2,2,4)
hold on
histogram(T_un_Numerical, 20, 'EdgeColor','b','FaceColor','b','normalization', 'pdf'); % numerical result
t_range = dt: dt: max(T_un_Numerical); 
plot(t_range, exp(- mu * t_range) * mu, 'r', 'linewidth', 2); % theoretical result; multiplying nu to make it to be a PDF
box on
set(gca,'fontsize',12)
title('Switching time from un to st states')
xlabel('t')
ylabel('Probability')