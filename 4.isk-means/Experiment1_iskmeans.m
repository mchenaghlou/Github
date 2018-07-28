clear all

% load dataset
load('./exp1_smiley_dataset_time_series.mat');


% call iskMeans algorithm
[ V, ClusterIndices, ~, ~, ~, ~, ~, ~, ~] = isKmeans(dataset(:, 1:2) , ...
    0, 0.99, 0.99, 0, './Experiment_1.avi');

% plot the dataset with prototypes.
figure('rend','painters','pos',[300 600 400 300])
hold on
for i = 1:max(ClusterIndices)
    plot(dataset(ClusterIndices==i, 1), dataset(ClusterIndices==i, 2), '.')
end

for i = 1:size(V, 1)
    plot(V(i, 1), V(i, 2), '*k', 'markersize', 12,'linewidth', 2)
end
 


axis equal
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

ylabel('Feature 2', 'FontSize',10,'FontWeight','bold')
xlabel('Feature 1', 'FontSize',10,'FontWeight','bold')

text(-10,50,'X_1', 'FontSize',12, 'FontWeight', 'bold')
text(14,50,'X_2', 'FontWeight', 'bold')
text(-35,7,'X_3', 'FontSize',12, 'FontWeight', 'bold')
text(-12,5,'X_4', 'FontSize',12, 'FontWeight', 'bold')
text(5,5,'X_5', 'FontSize',12, 'FontWeight', 'bold')
text(22,5,'X_6', 'FontSize',12, 'FontWeight', 'bold')
text(38,8,'X_7', 'FontSize',12, 'FontWeight', 'bold')

% calculate the Xie-Beni index value 
XB_batch(dataset, V)

% calculate the Adjusted Rand Index value 
rand_index(ClusterIndices, labels,'adjusted')

% calculate the Normalized Mutual Information value 
FindNMI(ClusterIndices', labels')