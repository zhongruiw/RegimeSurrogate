function [u_3R,h_W_3R,T_C_3R,T_E_3R,tau,I,s1,s2,s3,noise,n,sgm,p,xx] = reference_model(dt,T)
    %rng(22)
    rng(11)
    
    factor = 0.6; %check factor = 0.4, 0.5, 0.7 and 0.9
    
    [p,xx] = ksdensity([randn(1,15000)*0.6+0.5, randn(1,12500)*0.6+3],-0.5:0.1:4.5);
    xx = xx(6:end-5);
    p = p(6:end-5);
    p = p/trapz(xx,p);p=p*0+1/4;
    m = trapz(xx,xx.*p);
    
    lambda = 2/60;0.1;0.05;
    n = length(xx);
    Phi = zeros(1,n);
    sgm = zeros(1,n);
    for i = 2:n
        Phi(i) = trapz(xx(1:i), (xx(1:i)-m) .* p(1:i));
        sgm(i) = real(sqrt(2/p(i) * (-lambda * Phi(i))));
    end
    
    N = round(T/dt); % total numerical integration steps
    noise = randn(6,N);
    
    I = zeros(1,N);
    for i = 2:N
        temp = round(I(i-1)*10) + 1;
        if temp<0
            temp = 1;
        elseif temp > n
            temp = n;
        end
        sgm_x = sgm(temp);
        I(i) = I(i-1) + (-lambda * (I(i-1) - m) * dt) + sgm_x * noise(6,i-1) *sqrt(dt);
    end
    
    % Parameters of the deterministic part
    c = 1; gamma = 0.75 * factor; r = 0.25 * factor;  alpha_2 = 0.125 * factor; alpha_1 = alpha_2/2;  b_0 = 2.5;
    mu = 0.5; %sigma = 0.3;
    
    u_3R = zeros(1,N);
    h_W_3R = zeros(1,N);
    T_C_3R = zeros(1,N);
    T_E_3R = zeros(1,N);
    s1 = zeros(1,N);
    s2 = zeros(1,N);
    s3 = zeros(1,N);
    tau = zeros(1,N);
    d_tau = 2;
    u_3R(1) = 0;
    h_W_3R(1) = 0.0;
    
    
    Cu = 0.03  * factor;
    
    % Markov process parameters
    p_AtoB = 0.9;  % Probability of going from A to B
    p_BtoA = 0.9;  % Probability of going from B to A
    sigma_A = 0.1;  % Sigma value in state A
    sigma_B = 2.0;  % Sigma value in state B
    currentState = 'A'; % Initial state
    
    % Numerical integration
    for i = 2: N
        % sigma = I(i-1)/5 * factor;

        % Update state based on Markov process
        if currentState == 'A'
            if rand > p_AtoB
                currentState = 'B';
            end
        else  % currentState == 'B'
            if rand > p_BtoA
                currentState = 'A';
            end
        end

        % Set sigma based on current state
        if currentState == 'A'
            sigma = sigma_A;
        else  % currentState == 'B'
            sigma = sigma_B;
        end

        c1 = 26*(T_C_3R(i-1)+ 0.1).^2 + 0.95;  c1=c1*(1+0.4*sin(i*dt*2*pi/6-0*pi/6));
        c1 = c1  * factor;
        
        sigma_tau = 0.9* (tanh(4.5*T_C_3R(i-1)) + 1) *(1 + 0.25*1*cos(i*dt*2*pi/6+0*2*pi/6));
        beta_E = 1.0+(1-I(i-1)/5);beta_E=beta_E/10*1.6 * sqrt(factor);
        beta_u = -0.2*beta_E;
        beta_h = -0.4*beta_E;
        beta_C = 0.8*beta_E;
        
        c2= 1.5 * factor*( 1 + 0.4*sin(i*dt*2*pi/6 + 2*pi/6) + 0.2* sin(2*i*dt*2*pi/6 + 2*pi/6));

        
        u_3R(i) = u_3R(i-1) + ( - r * u_3R(i-1) - alpha_1 * b_0 * mu / 2 * (T_C_3R(i-1) + T_E_3R(i-1)) ) * dt + 0.04 * sqrt(factor) * sqrt(dt) * noise(1,i-1) + beta_u * tau(i-1) * dt;
        h_W_3R(i) = h_W_3R(i-1) + ( - r * h_W_3R(i-1) - alpha_2 * b_0 * mu / 2 * (T_C_3R(i-1) + T_E_3R(i-1)) ) * dt + 0.02 * sqrt(factor) * sqrt(dt) * noise(2,i-1) + beta_h * tau(i-1) * dt;
        T_C_3R(i) = T_C_3R(i-1) + ( (gamma * b_0 * mu / 2 - c1) * T_C_3R(i-1) + gamma * b_0 * mu / 2 * T_E_3R(i-1) + gamma * h_W_3R(i-1) + sigma * u_3R(i-1) + Cu ) * dt + beta_C * tau(i-1) * dt+ 0.02*2 * sqrt(factor) * sqrt(dt) * noise(3,i-1);
        T_E_3R(i) = T_E_3R(i-1) + ( gamma * h_W_3R(i-1) + (1 * gamma * b_0 * mu / 2 +gamma * b_0 * mu - c2) * T_E_3R(i-1) + (1 * gamma * b_0 * mu / 2 - gamma * b_0 * mu) * T_C_3R(i-1)) * dt + beta_E * tau(i-1) * dt+ 0.03 *1* sqrt(dt)* sqrt(factor) * noise(4,i-1);
        tau(i) = tau(i-1) - d_tau * tau(i-1) * dt + sigma_tau * noise(5,i-1) * sqrt(dt);
    
        s1(i-1) = sin(pi*i*dt/3-0*pi/6);
        s2(i-1) = sin(pi*i*dt/3 + 2*pi/6);
        s3(i-1) = sin(2*pi*i*dt/3 + 2*pi/6);
    end
end