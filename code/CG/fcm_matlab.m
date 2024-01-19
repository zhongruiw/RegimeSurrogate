clear all

rng(2023)

c=16;
% save('cluster_idx.mat', "cluster_idx")

mydata = load('../../data/naturerun_pwell_step.mat');
your_data =  mydata.data';
[data, N] = size(mydata.data);  % Assuming your_data is your 3xN matrix

% Initialize cluster centers 
% centers = your_data(randperm(N, c), :); % randomly
% [idx,centers] = kmeans(your_data,c,'MaxIter',500); % by kmeans
centers = zeros(16,3);
centers(1,:) = [1.8,1.8,1.8];
centers(2,:) = [-1.8,1.8,1.8];
centers(3,:) = [1.8,-1.8,1.8];
centers(4,:) = [1.8,1.8,-1.8];
centers(5,:) = [-1.8,-1.8,1.8];
centers(6,:) = [1.8,-1.8,-1.8];
centers(7,:) = [-1.8,1.8,-1.8];
centers(8,:) = [-1.8,-1.8,-1.8];
centers(9,:) = [1.0,1.0,1.0];
centers(10,:) = [-1.0,1.0,1.0];
centers(11,:) = [1.0,-1.0,1.0];
centers(12,:) = [1.0,1.0,-1.0];
centers(13,:) = [-1.0,-1.0,1.0];
centers(14,:) = [1.0,-1.0,-1.0];
centers(15,:) = [-1.0,1.0,-1.0];
centers(16,:) = [-1.0,-1.0,-1.0];


% Perform FCM clustering
options = fcmOptions(...
    NumClusters=16,...
    Exponent=2,...
    ClusterCenters=centers,...
    Verbose=false);
[center, U] = fcm(your_data, options);

[max_val, cluster_idx] = max(U);

cluster_idx = cluster_idx';
% Create a color map for the clusters
% colors = jet(c);
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
    cluster_indices = find(cluster_idx == i);
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
title('Fuzzy C-Means Clustering of Data');
grid on;
hold off;

figure;
scatter3(center(:,1),center(:,2),center(:,3),80)
xlabel('x');
ylabel('y');
zlabel('z');
title('Centers of FCM');

figure;
scatter3(centers(:,1),centers(:,2),centers(:,3),80)
xlabel('x');
ylabel('y');
zlabel('z');
title('initial centers');