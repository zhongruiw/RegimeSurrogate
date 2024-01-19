% clc 
% clear all
% 
% load reg_te_tc_hw_ssta
% load lon_sst
% load obs_data.mat
% 
% % loaded_data = load('/Users/ree/Documents/Research/RegimeSurrogate/code/ClusteringDynamical/data/potential_well_step_Ludo.mat');
% 
% dt = 0.002; % time step
% T = 12000; % total time length, i.e., 1200*2=2400 months=200 years
% N = round(T/dt); % total numerical integration steps
% 
% [u_3R,h_W_3R,T_C_3R,T_E_3R,tau,I,s1,s2,s3,noise,n,sgm,p,xx] = reference_model(dt,T);
% 
% u_lg = 1.5;
% h_lg = 150;
% t_lg = 7.5;
% tau_lg = 5;
% 
% % upsampling 
% res = 50; % about 6 days
% T_C = T_C_3R(1:res:end);
% T_E = T_E_3R(1:res:end);
% h_W = h_W_3R(1:res:end);
% u_ = u_3R(1:res:end);
% tau_ = tau(1:res:end);
% T_C=T_C-mean(T_C);
% T_E=T_E-mean(T_E);
% h_W=h_W-mean(h_W);
% u_=u_-mean(u_);
% tau_=tau_-mean(tau_);
% T_C = T_C * 7.5;
% T_E = T_E * 7.5;
% h_W = h_W * 150;
% u_ = u_ * 1.5;
% tau_ = tau_ * tau_lg;
% h_W_obs_new = h_W_obs(1:end-12)';
% u_obs_new = u_obs';
% T_E_obs_new = T_E_obs(30*12+1:end-12)';
% T_C_obs_new = T_C_obs(30*12+1:end-12)';
% 
% k_dt=0.5/dt;
% nn=1;
% T_C_3R=smooth(T_C_3R,nn*k_dt);
% T_C_3R=T_C_3R(1:k_dt:end);
% T_E_3R=smooth(T_E_3R,nn*k_dt);
% T_E_3R=T_E_3R(1:k_dt:end);
% h_W_3R=smooth(h_W_3R,nn*k_dt);
% h_W_3R=h_W_3R(1:k_dt:end);
% u_3R=smooth(u_3R,nn*k_dt);
% u_3R=u_3R(1:k_dt:end);
% I=smooth(I,nn*k_dt);
% I=I(1:k_dt:end);
% tau_3R=smooth(tau,nn*k_dt);
% tau_3R=tau_3R(1:k_dt:end);
% T_C_3R=T_C_3R-mean(T_C_3R);
% T_E_3R=T_E_3R-mean(T_E_3R);
% h_W_3R=h_W_3R-mean(h_W_3R);
% u_3R=u_3R-mean(u_3R);
% tau_3R=tau_3R-mean(tau_3R);
% T_C_3R = T_C_3R * 7.5;
% T_E_3R = T_E_3R * 7.5;
% h_W_3R = h_W_3R * 150;
% u_3R = u_3R * 1.5;
% tau_3R = tau_3R * tau_lg;
% 
% % write_data('enso_sigma2state','T_C',T_C,'T_E',T_E,'h_W',h_W,'u',u_,'tau',tau_,'Hov',hov,'dt',dt*res,'res',res,'T_C_obs',T_C_obs_new','T_E_obs',T_E_obs_new','h_W_obs',h_W_obs_new','u_obs',u_obs_new');
% 
% colors = [
%     0 0 0;    % Black
%     1 0 1;    % Magenta
%     0 0 1;    % Blue
%     1 0.5 0;  % Orange
%     0.8 0.8 0;  % Olive
%     1 0 0;    % Red
%     0.5 0.5 0.5;  % Gray
%     0 0.8 0.8;  % Teal
%     0 1 0;    % Green
%     1 1 0;    % Yellow
%     0.8 0 0.8;  % Pink
%     0 1 1;    % Cyan
%     0.5 0 1;  % Purple
%     0 0.5 1;  % Light Blue
%     0.5 1 0;  % Light Green
%     0 0.3 0;  % Dark Green
% ];
% 
% % select features
% data = [T_C_3R, T_E_3R];%, h_W_3R, u_3R,tau_3R];  % need tune here-------------------------
% [K,nfeature] = size(data);
% % Xtest = [T_C_obs_new, T_E_obs_new]; % need tune here---------------------
% % [K_obs,nfeature] = size(Xtest);
% 
% % standardize data
% data = zscore(data);
% % Xtest = zscore(Xtest);
% 
% % clustering
% cluster_method = 'Ludos';
% if strcmp(cluster_method, 'Kmeans')
%     c = [2,3,4];
% elseif strcmp(cluster_method, 'Ludos')
%     % load Ludo's clustering results
%     loaded_data = load('/Users/ree/Documents/Research/RegimeSurrogate/code/ClusteringDynamical/data/enso_sigma2state_TcTe_Ludo_standard.mat');
%     % Access your variables
%     idxs = loaded_data.X_LN_array;
%     % for i = 1:size(idxs, 1)
%     %     idxs{i} = idxs{i}(:, 2);
%     % end
%     ncs = cellfun(@max,idxs);
%     tscales = loaded_data.tau;
%     c = ncs(1:8); % choose ncenters
% end
% 
% finecluster = load('/Users/ree/Documents/Research/RegimeSurrogate/code/ClusteringDynamical/data/enso_TcTe_finecluster.mat');
% X = finecluster.X;
% Xc = finecluster.Xc;

tmp = idxs{4,1}(X==1222);
% tmp(1)

% finecluster = load('/Users/ree/Documents/Research/RegimeSurrogate/code/ClusteringDynamical/data/potential_well_finecluster.mat');
% X = finecluster.X;
% Xc = finecluster.Xc;
% 
% % tmp = idxs{4,1}(X==6);
% % tmp(1)
