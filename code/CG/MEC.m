clear all
rng(2023)

mydata = load('../../data/naturerun_pwell.mat');
data =  mydata.data';
K = 16;
[idx,C] = kmeans(data,K,'MaxIter',500);
% cluster_idx = load("cluster_idx.mat");
% label_old = cluster_idx.cluster_idx;
gama = 0.3;
[label,iter_FCM, para_miu, NegativeLogLikelihood, responsivity]=MEC_kailugaji(data, K, C, gama);

% save data
mec_result = struct();
mec_result.label = label;  % Replace 'label' with your actual variable
mec_result.iter_FCM = iter_FCM;  % Replace 'iter_FCM' with your actual variable
mec_result.para_miu = para_miu;  % Replace 'para_miu' with your actual variable
mec_result.NegativeLogLikelihood = NegativeLogLikelihood;  % Replace 'NegativeLogLikelihood' with your actual variable
mec_result.responsivity = responsivity;  % Replace 'responsivity' with your actual variable
save('../../data/mec_pwell_step.mat', 'mec_result');

c =K;
your_data = data;
colors = [
    0 0 1;    % Blue
    0 1 0;    % Green
    1 0 0;    % Red
    0 1 1;    % Cyan
    1 0 1;    % Magenta
    1 1 0;    % Yellow
    0 0.5 1;  % Light Blue
    0.5 0 1;  % Purple
    0.5 1 0;  % Light Green
    1 0.5 0;  % Orange
    0 0.8 0.8;  % Teal
    0.8 0 0.8;  % Pink
    0.8 0.8 0;  % Olive
    0 0.3 0;  % Dark Green
    0.5 0.5 0.5;  % Gray
    0 0 0;    % Black
];

% Create a figure for the scatter plot
figure;
hold on;

% Plot data points for each cluster in a different color
for i = 1:c
    cluster_indices = find(label == i);
    cluster_data = your_data(cluster_indices,:);
    x = cluster_data(:,1);
    y = cluster_data(:,2);
    z = cluster_data(:,3);
    scatter3(x, y, z, 1, colors(i, :), 'filled');
end

% Customize the plot
xlabel('x');
ylabel('y');
zlabel('z');
title('Maximum Entropy Clustering of Data');
grid on;
hold off;

figure;
scatter3(C(:,1),C(:,2),C(:,3),80)
xlabel('x');
ylabel('y');
zlabel('z');
title('Center of Kmeans');

function [label,iter_FCM, para_miu, NegativeLogLikelihood, responsivity]=MEC_kailugaji(data, K, C, gama)
% Input:
% K: number of cluster
% data: dataset, N*D
% label_old: initializing label. N*1
% Output:
% label: results of cluster. N*1
% iter_FCM: iterations
% Written by kailugaji. (wangrongrong1996@126.com)
format long 

%% initializing parameters
esp=1e-6;  % stopping criterion for iteration
max_iter=100;  % maximum number of iterations 
fitness=zeros(max_iter,1);
[data_num,data_dim]=size(data);
distant=zeros(data_num, K);
responsivity=zeros(data_num,K);
para_miu=zeros(K, data_dim);
%% initializing the cluster center
% % random
% rand_array=randperm(data_num);   
% para_miu=data(rand_array(1:K), :); 
para_miu = C;
% for k=1:K
%     X_k=data(label_old==k, :); 
%     para_miu(k, :)=mean(X_k); % the center of each cluster  
% end
%% Maximum entropy clustering algorithm
for t=1:max_iter
    % (X-para_miu)^2=X^2+para_miu^2-2*para_miu*X'. data_num*K
    for k=1:K
        distant(:,k)=sum((data-repmat(para_miu(k, :), data_num, 1)).^2,2); %N*1
    end
    % update membership. data_num*K
    R_up=exp(-distant./gama);  
    responsivity= R_up./repmat(sum(R_up,2), 1, K);
    % update center. K*data_dim
    miu_up=(responsivity')*data;  
    para_miu=miu_up./((sum(responsivity))'*ones(1,data_dim));
    % object function
    fitness(t)=sum(sum(responsivity.*distant))+gama.*sum(sum((responsivity.*log(responsivity+eps))));
    if t>1  
        if abs(fitness(t)-fitness(t-1))<esp
            break;
        end
    end
end
iter_FCM=t;  % iterations
NegativeLogLikelihood=fitness(iter_FCM);
%% clustering
[~,label]=max(responsivity,[],2);

end
