clc 
clear all

load reg_te_tc_hw_ssta
load lon_sst
load obs_data.mat

% loaded_data = load('/Users/ree/Documents/Research/RegimeSurrogate/code/ClusteringDynamical/data/potential_well_step_Ludo.mat');

dt = 0.002; % time step
T = 12000; % total time length, i.e., 1200*2=2400 months=200 years
N = round(T/dt); % total numerical integration steps

[u_3R,h_W_3R,T_C_3R,T_E_3R,tau,I,s1,s2,s3,noise,n,sgm,p,xx] = reference_model(dt,T);

write_data('../data/enso_sigma2state_referrun1','T_C',T_C_3R,'T_E',T_E_3R,'h_W',h_W_3R,'u',u_3R,'tau',tau,'dt',dt,'I',I,'s1',s1,'s2',s2,'s3',s3,'noise',noise,'n',n,'sgm',sgm,'p',p,'xx',xx);

u_lg = 1.5;
h_lg = 150;
t_lg = 7.5;
tau_lg = 5;

% upsampling 
res = 50; % about 6 days
T_C = T_C_3R(1:res:end);
T_E = T_E_3R(1:res:end);
h_W = h_W_3R(1:res:end);
u_ = u_3R(1:res:end);
tau_ = tau(1:res:end);
T_C=T_C-mean(T_C);
T_E=T_E-mean(T_E);
h_W=h_W-mean(h_W);
u_=u_-mean(u_);
tau_=tau_-mean(tau_);
T_C = T_C * 7.5;
T_E = T_E * 7.5;
h_W = h_W * 150;
u_ = u_ * 1.5;
tau_ = tau_ * tau_lg;
h_W_obs_new = h_W_obs(1:end-12)';
u_obs_new = u_obs';
T_E_obs_new = T_E_obs(30*12+1:end-12)';
T_C_obs_new = T_C_obs(30*12+1:end-12)';

k_dt=0.5/dt;
nn=1;
T_C_3R=smooth(T_C_3R,nn*k_dt);
T_C_3R=T_C_3R(1:k_dt:end);
T_E_3R=smooth(T_E_3R,nn*k_dt);
T_E_3R=T_E_3R(1:k_dt:end);
h_W_3R=smooth(h_W_3R,nn*k_dt);
h_W_3R=h_W_3R(1:k_dt:end);
u_3R=smooth(u_3R,nn*k_dt);
u_3R=u_3R(1:k_dt:end);
I=smooth(I,nn*k_dt);
I=I(1:k_dt:end);
tau_3R=smooth(tau,nn*k_dt);
tau_3R=tau_3R(1:k_dt:end);
T_C_3R=T_C_3R-mean(T_C_3R);
T_E_3R=T_E_3R-mean(T_E_3R);
h_W_3R=h_W_3R-mean(h_W_3R);
u_3R=u_3R-mean(u_3R);
tau_3R=tau_3R-mean(tau_3R);
T_C_3R = T_C_3R * 7.5;
T_E_3R = T_E_3R * 7.5;
h_W_3R = h_W_3R * 150;
u_3R = u_3R * 1.5;
tau_3R = tau_3R * tau_lg;

res_lon = 10;
hov = zeros(length(lon_sst)/res_lon, length(T_C));
for i =1:length(T_C)
    hov_temp = Hovmoller(lon_sst, T_E(i), T_C(i),reg_te_tc_hw_ssta);
    hov(:,i) = hov_temp(1:res_lon:end);
end

% write_data('enso_sigma2state','T_C',T_C,'T_E',T_E,'h_W',h_W,'u',u_,'tau',tau_,'Hov',hov,'dt',dt*res,'res',res,'T_C_obs',T_C_obs_new','T_E_obs',T_E_obs_new','h_W_obs',h_W_obs_new','u_obs',u_obs_new');

colors = [
    0 0 0;    % Black
    1 0 1;    % Magenta
    0 0 1;    % Blue
    1 0.5 0;  % Orange
    0.8 0.8 0;  % Olive
    1 0 0;    % Red
    0.5 0.5 0.5;  % Gray
    0 0.8 0.8;  % Teal
    0 1 0;    % Green
    1 1 0;    % Yellow
    0.8 0 0.8;  % Pink
    0 1 1;    % Cyan
    0.5 0 1;  % Purple
    0 0.5 1;  % Light Blue
    0.5 1 0;  % Light Green
    0 0.3 0;  % Dark Green
];

% select features
data = [T_C_3R, T_E_3R];%, h_W_3R, u_3R,tau_3R];  % need tune here-------------------------
[K,nfeature] = size(data);
% Xtest = [T_C_obs_new, T_E_obs_new]; % need tune here---------------------
% [K_obs,nfeature] = size(Xtest);

% standardize data
data = zscore(data);
% Xtest = zscore(Xtest);

% clustering
cluster_method = 'Kmeans';
if strcmp(cluster_method, 'Kmeans')
    c = [2,3,4];
elseif strcmp(cluster_method, 'Ludos')
    % load Ludo's clustering results
    loaded_data = load('/Users/ree/Documents/Research/RegimeSurrogate/code/ClusteringDynamical/data/enso_sigma2state_TcTe_Ludo_standard_2cluster.mat');
    % Access your variables
    idxs = loaded_data.X_LN_array;
    for i = 1:size(idxs, 1)
        idxs{i} = idxs{i}(:, 2);
    end
    ncs = cellfun(@max,idxs);
    tscales = loaded_data.tau;
    c = ncs(1:8); % choose ncenters
end

clusterIndices = cell(length(c), 1); % store the indices of each cluster
centroids = cell(length(c), 1); % store the center of each cluster

% clustering based on boreal_winter or not
boreal_winter = false;

% zoom in to see each month's label or not
zoomin = false;

% Loop over the number of clusters and perform k-means clustering
for nc = 1:length(c)
    if boreal_winter
        boreal_winter_means = calculateBorealWinterMeans(data);
        % boreal_winter_means_test = calculateBorealWinterMeans(Xtest);

        % Perform k-means clustering
        idx = zeros(K,1);
        [idx_bw, C] = kmeans(boreal_winter_means, c(nc));
        for year = 1:K/12-1
            idx(year*12:year*12+2) = idx_bw(year);
        end

        % % test on observations
        % idx_test = zeros(K_obs,1);
        % [~,idx_bw_test] = pdist2(C,boreal_winter_means_test,'euclidean','Smallest',1);
        % for year = 1:K_obs/12-1
        %     idx_test(year*12:year*12+2) = idx_bw_test(year);
        % end
    else
        if strcmp(cluster_method, 'Kmeans')
            % Perform k-means clustering
            [idx, C] = kmeans(data, c(nc));
        elseif strcmp(cluster_method, 'Ludos')
            repeated_idx = repmat(idxs{nc,1}, 1,res);
            idx = reshape(repeated_idx',[],1);
            idx = idx(1:k_dt:end);
        end

        % % test on observations --------------------------- may need tune here
        % [~,idx_test] = pdist2(C,Xtest,'euclidean','Smallest',1);
    end

    % centers{nc} = C;
    idxClusters{nc} = idx;

    %% Model_hov_particular_events
    if zoomin
        year_add = 2;
        t_temp = 1980+1/12:1/12:1980+year_add;%20 years
        window = 1;
        total_loop = (year_add) * 12/window-1;
    else
        year_add = 20;
        t_temp = 1980+1/12:1/12:1980+year_add;%20 years
        window = 12;
        total_loop = (year_add-1) * 12/window;
    end
    LL = length(t_temp);

    fig = figure('color','white');
    set(gcf,'unit','centimeters','position',[10 5 18 15])
    colormap(jet)
    for j = 1:5
        if j == 1
            range_model = 1+12*220:LL+12*220;
    %                 range_model = 1+12*900:LL+12*900;
        elseif j == 2
            range_model = 1+12*240:LL+12*240;
    %                 range_model = 1+12*920:LL+12*920;
        elseif j == 3
            range_model = 1+12*1410:LL+12*1410;
    %                 range_model = 1+12*940:LL+12*940;
        elseif j == 4
            range_model = 1+12*750:LL+12*750;
    %                 range_model = 1+12*960:LL+12*960;
        elseif j == 5
            range_model = 1+12*900:LL+12*900;
        end
        range_tau = range_model(1)*k_dt:range_model(end)*k_dt;
        t_model = range_model/12;
        tau_model = range_tau/k_dt/12;
        Hov=zeros(length(lon_sst),length(range_model));

        xx=[ones(size(T_E_3R(range_model))),T_E_3R(range_model),T_C_3R(range_model)];
        for i=1:length(lon_sst)
            Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
        end

        subplot(1,5,j)  % may tune here -----------------------------------
        [xx,yy] = meshgrid(t_model,lon_sst);
        contourf(yy,xx,Hov,30,'linestyle','none')
        hold on
        plot([180 180],[t_model(1) t_model(end)],'m--','linewidth',2);
        temp_tau = range_tau(1:k_dt/3:end);
        plot(180+tau(temp_tau)*20,tau_model(1:k_dt/3:end),'k','linewidth',0.5);
        for k = 1:total_loop
            for ilabel = 1:c(nc)+1
                if idx(range_model(1)+k*window) == ilabel
                % if all(idx(range_model(1)+k*window-1:range_model(1)+k*window+1) == ilabel)
                    if zoomin
                        plot([120,120],[t_model(window*k)-0.02,t_model(window*k)+0.02],'color',colors(ilabel,:),'linewidth',10)
                    else
                        plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'color',colors(ilabel,:),'linewidth',10)
                    end
                end
            end

        end
        %     colorbar
        set(gca,'xlim',[120 280]);
        set(gca,'xtick',120:60:280);
        % if j == 3
        %     set(gca,'yticklabel',402:2:420);
        % end
        % if j == 4
        %     set(gca,'yticklabel',442:2:460);
        % end
        % if j == 5
        %     set(gca,'yticklabel',502:2:520);
        % end

        xlabel('longitude');
        set(gca,'fontsize',9,'linewidth',1);
        %     title('Evolutions of SSTA and wind burst','fontsize',9);
        caxis([-3,3])
        %     colorbar('eastoutside','position',[.92,0.12,0.005,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');
        if j == 1
            ylabel('model year');
            %         text(460,241,'Hovmoller diagram of the standard run ','fontsize',16)
        end
    end

    % % test on observations ------------------------------ may tune here
    % j = 6;
    % if zoomin
    %     t_obs = 1999+1/12:1/12:1999+year_add;%20 years
    % else
    %     t_obs = 1999+1/12:1/12:1999+year_add;%20 years
    % end
    % LL = length(t_obs);
    % range_obs = 1+12*(1999-1980):LL+12*(1999-1980);
    % range_tau = range_obs(1)*k_dt:range_obs(end)*k_dt;
    % tau_obs = range_tau/k_dt/12;
    % Hov=zeros(length(lon_sst),length(range_obs));
    % 
    % xx=[ones(size(T_E_obs_new(range_obs))),T_E_obs_new(range_obs),T_C_obs_new(range_obs)];
    % for i=1:length(lon_sst)
    %     Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
    % end
    % subplot(1,6,j)
    % [xx,yy] = meshgrid(t_obs,lon_sst);
    % contourf(yy,xx,Hov,30,'linestyle','none')
    % hold on
    % plot([180 180],[t_obs(1) t_obs(end)],'m--','linewidth',2);
    % temp_tau = range_tau(1:k_dt/3:end);
    % % plot(180+tau(temp_tau)*20,tau_obs(1:k_dt/3:end),'k','linewidth',0.5);
    % for k = 1:total_loop
    %     for ilabel = 1:nc+1
    %         if idx_test(range_obs(1)+k*window) == ilabel
    %         % if all(idx_test(range_obs(1)+k*window-1:range_obs(1)+k*window+1) == ilabel)
    %             plot([120,120],[t_obs(1-4+window*k),t_obs(1+1+window*k)],'color',colors(ilabel,:),'linewidth',10)
    %         end
    %     end
    % end
    % %     colorbar
    % set(gca,'xlim',[120 280]);
    % set(gca,'xtick',120:60:280);
    % xlabel('longitude');
    % set(gca,'fontsize',9,'linewidth',1);
    % %     title('Evolutions of SSTA and wind burst','fontsize',9);
    % caxis([-3,3])

    % sgtitle(sprintf('Hovmoller diagrams of %d centers clustering based on T_C,T_E, tscale=%.2f', c(nc),tscales(nc)), 'FontSize', 12); % need tune here-------------------------
    sgtitle(sprintf('Hovmoller diagrams of %d centers clustering based on T_C,T_E', c(nc)), 'FontSize', 12); % need tune here-------------------------
    colorbar('eastoutside','position',[.92,0.12,0.01,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');

    % save figure
    filename = sprintf('./figure/Hov_%d_centers_TcTe_standard_sig2state.png', c(nc)); % need tune here-------------------------
    % filename = sprintf('./figure/Hov_%d_centers_TcTe_standard_sig2state_tscale%.3f_2cluster.png', c(nc),tscales(nc)); % need tune here-------------------------
    print(fig, filename, '-dpng', '-r150');

end

% % Plot each ENSO events evolution
% fig = figure('color','white');
% set(gcf,'unit','centimeters','position',[5 15 12 10])
% hold on; % Keep the same plot to overlay new lines
% colors = lines(5); % Generate a set of colors for the lines
% 
% year_add = 4;
% t_temp = 1980+1/12:1/12:1980+year_add;%20 years
% LL = length(t_temp);
% window = 12;
% total_loop = (year_add-1) * 12/window;
% 
% for j = 1:5
%     range_model = 1+12*(750+(j-1)*year_add):LL+12*(750+(j-1)*year_add);
% 
%     t_model = range_model/12;
% 
%     % Generate some example data for each event
%     Te_data = T_E_3R(range_model);
%     Tc_data = T_C_3R(range_model);
% 
%     % Plot the event data
%     plot(Te_data, Tc_data, 'o-', 'Color', colors(j,:), 'LineWidth', 1.0, ...
%          'MarkerFaceColor', colors(j,:), 'DisplayName', sprintf('%3.0f-%d', t_model(1), t_model(end)));
% end
% 
% % Add a legend
% legend('show');
% 
% % Labels and title
% xlabel('T_E (°C)');
% ylabel('T_C (°C)');
% title('Evolution of ENSO events');
% 
% % Draw zero lines
% line(xlim, [0 0], 'Color', 'black', 'LineWidth', 0.5,'HandleVisibility','off'); % Horizontal line
% line([0 0], ylim, 'Color', 'black', 'LineWidth', 0.5,'HandleVisibility','off'); % Vertical line
% 
% % Grid and hold off
% grid on;
% hold off;
% filename = sprintf('./figure/TcTe_phase_sig2.png'); % need tune here-------------------------
% print(fig, filename, '-dpng', '-r150');

% % manually defined ENSO types
% year_add = 20;
% t_temp = 1980+1/12:1/12:1980+year_add;%20 years
% LL = length(t_temp);
% 
% window = 12;
% total_loop = (year_add-1) * 12/window;
% fig = figure('color','white');
% set(gcf,'unit','centimeters','position',[10 5 18 15])
% colormap(jet)
% for j = 1:5
%     if j == 1
%         range_model = 1+12*220:LL+12*220;
% %                 range_model = 1+12*900:LL+12*900;
%     elseif j == 2
%         range_model = 1+12*240:LL+12*240;
% %                 range_model = 1+12*920:LL+12*920;
%     elseif j == 3
%         range_model = 1+12*1410:LL+12*1410;
% %                 range_model = 1+12*940:LL+12*940;
%     elseif j == 4
%         range_model = 1+12*750:LL+12*750;
% %                 range_model = 1+12*960:LL+12*960;
%     elseif j == 5
%         range_model = 1+12*900:LL+12*900;
%     end
%     range_tau = range_model(1)*k_dt:range_model(end)*k_dt;
%     t_model = range_model/12;
%     tau_model = range_tau/k_dt/12;
%     Hov=zeros(length(lon_sst),length(range_model));
% 
%     xx=[ones(size(T_E_3R(range_model))),T_E_3R(range_model),T_C_3R(range_model)];
%     for i=1:length(lon_sst)
%         Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
%     end
%     subplot(1,6,j)
%     [xx,yy] = meshgrid(t_model,lon_sst);
%     contourf(yy,xx,Hov,30,'linestyle','none')
%     hold on
%     plot([180 180],[t_model(1) t_model(end)],'m--','linewidth',2);
%     temp_tau = range_tau(1:k_dt/3:end);
%     plot(180+tau(temp_tau)*20,tau_model(1:k_dt/3:end),'k','linewidth',0.5);
%     for k = 1:total_loop
%         if mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 1.0 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
%             plot([120,120],[t_model(1-4+k*window),t_model(1+1+window*k)],'k','linewidth',10)
%         elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
%             plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'b','linewidth',10)
%         elseif mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
%             plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'m','linewidth',10)
%         elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5 || mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5
%             plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'color',[255 97 0]/255,'linewidth',10)
%         end
%     end
%     %     colorbar
%     set(gca,'xlim',[120 280]);
%     set(gca,'xtick',120:60:280);
% %     if j == 3
% %         set(gca,'yticklabel',402:2:420);
% %     end
% %     if j == 4
% %         set(gca,'yticklabel',442:2:460);
% %     end
% %     if j == 5
% %         set(gca,'yticklabel',502:2:520);
% %     end
% 
%     xlabel('longitude');
%     set(gca,'fontsize',9,'linewidth',1);
%     %     title('Evolutions of SSTA and wind burst','fontsize',9);
%     caxis([-3,3])
%     %     colorbar('eastoutside','position',[.92,0.12,0.005,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');
%     if j == 1
%         ylabel('model year');
%         %         text(460,241,'Hovmoller diagram of the standard run ','fontsize',16)
%     end
% end
% 
% % test on observations
% j = 6;
% t_obs = 1999+1/12:1/12:2019;
% LL = length(t_obs);
% range_obs = 1+12*(1999-1980):LL+12*(1999-1980);
% range_tau = range_obs(1)*k_dt:range_obs(end)*k_dt;
% tau_obs = range_tau/k_dt/12;
% Hov=zeros(length(lon_sst),length(range_obs));
% 
% xx=[ones(size(T_E_obs_new(range_obs))),T_E_obs_new(range_obs),T_C_obs_new(range_obs)];
% for i=1:length(lon_sst)
%     Hov(i,:)=xx*reg_te_tc_hw_ssta(i,1:3)';
% end
% subplot(1,6,j)
% [xx,yy] = meshgrid(t_obs,lon_sst);
% contourf(yy,xx,Hov,30,'linestyle','none')
% hold on
% plot([180 180],[t_obs(1) t_obs(end)],'m--','linewidth',2);
% temp_tau = range_tau(1:k_dt/3:end);
% % plot(180+tau(temp_tau)*20,tau_obs(1:k_dt/3:end),'k','linewidth',0.5);
% for k = 1:total_loop
%     if mean(T_E_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window)) > 1.0 && mean(T_E_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window))>mean(T_C_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window))
%         plot([120,120],[t_obs(1-4+k*window),t_obs(1+1+window*k)],'k','linewidth',10)
%     elseif mean(T_E_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window)) > 0.5 && mean(T_E_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window))>mean(T_C_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window))
%         plot([120,120],[t_obs(1-4+window*k),t_obs(1+1+window*k)],'b','linewidth',10)
%     elseif mean(T_C_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window)) > 0.5 && mean(T_C_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window))>mean(T_E_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window))
%         plot([120,120],[t_obs(1-4+window*k),t_obs(1+1+window*k)],'m','linewidth',10)
%     elseif mean(T_E_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window)) < -0.5 || mean(T_C_obs_new(range_obs(1)-1+k*window:range_obs(1)+1+k*window)) < -0.5
%         plot([120,120],[t_obs(1-4+window*k),t_obs(1+1+window*k)],'color',[255 97 0]/255,'linewidth',10)
%     end
% end
% %     colorbar
% set(gca,'xlim',[120 280]);
% set(gca,'xtick',120:60:280);
% xlabel('longitude');
% set(gca,'fontsize',9,'linewidth',1);
% %     title('Evolutions of SSTA and wind burst','fontsize',9);
% caxis([-3,3])
% 
% sgtitle('Hovmoller diagrams of manually defined ENSO','fontsize',12);
% colorbar('eastoutside','position',[.92,0.12,0.01,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');
% 
% % save figure
% filename = sprintf('./figure/Hov_manenso.png');
% print(fig, filename, '-dpng', '-r150');


% 

function boreal_winter_means = calculateBorealWinterMeans(data)
    [K,nfeature] = size(data);
    num_years = K / 12;
        
    % Preallocate array for boreal winter means
    boreal_winter_means = zeros(num_years-1, nfeature);
    
    % Loop through each year
    for year = 1:num_years
        % Indices for Dec-Jan-Feb
        dec_idx = (year - 1) * 12 + 12; % December of previous year
        jan_idx = (year - 1) * 12 + 13; % January of the current year
        feb_idx = (year - 1) * 12 + 14; % February of the current year
        
        % Check for out of bounds (for the last year)
        if jan_idx > K
            break;
        end

        % Extract boreal winter data
        boreal_winter_data = data([dec_idx, jan_idx, feb_idx],:);
        
        % Calculate mean of boreal winter data
        boreal_winter_means(year,:) = mean(boreal_winter_data,1);
    end
end

function write_data(filename, varargin)
    % Define the full path to the HDF5 file
    hdf5file = fullfile('data', [filename '.hdf5']);
    
    % Iterate over input arguments
    for k = 1:2:length(varargin)
        dataName = varargin{k}; % Name of the variable
        dataValue = varargin{k+1}; % Value of the variable
        
        % Create the dataset
        h5create(hdf5file, ['/' dataName], size(dataValue));
        
        % Write the data to the dataset
        h5write(hdf5file, ['/' dataName], dataValue);
    end
end

function hov = Hovmoller(lon_sst,T_e,T_c,reg_te_tc_hw_ssta)
    hov=zeros(length(lon_sst),1);
    xx=[1,T_e,T_c]; % xx should be (1,3) size
    for i=1:length(lon_sst)
        hov(i)=xx*reg_te_tc_hw_ssta(i,1:3)';
    end
end
