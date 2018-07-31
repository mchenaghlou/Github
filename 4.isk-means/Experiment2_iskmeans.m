clear all
% Dataset_smiley
% save('exp2_smiley_dataset.mat', 'dataset','labels')
load('exp2_smiley_dataset.mat')

% plot(dataset(:, 1), dataset(:, 2), '.');
tic
[ V, ClusterIndices, iCVs_ff, iCVs1_ff, added, removed ] = isKmeans(dataset , 2, 0.99, 0.99, 0, ...
    './Experiment_2.avi')
time = toc

figure('rend','painters','pos',[300 600 400 300])
hold on
for i = 1:max(ClusterIndices)
    plot(dataset(ClusterIndices==i, 1), dataset(ClusterIndices==i, 2), '.')
end
plot(V(:, 1), V(:, 2), '*k', 'markersize', 5, 'LineWidth',5)
axis equal

% xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

% yt = get(gca, 'YTick');
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



XB_batch(dataset, V);

[ClusterIndices, labels] = ReconstructLabels(ClusterIndices, labels, false);
rand_index(ClusterIndices, labels,'adjusted')

FindNMI(ClusterIndices', labels')

