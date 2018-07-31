clear all
filename = './GSA_dataset.mat';
load(filename);
dim = size(reduced_subset,2) - 1;
labels = reduced_subset(:, end);
dataset = reduced_subset(:, 1:end-1);

tic
% [1:15000,16007:17007,24008:25008,34009:40000]
[ V, ClusterIndices, iCVs_ff, iCVs1_ff, added, removed] = isKmeans(dataset, 0, 0.99, 0.99, 0, ...
    './4.isk-means/Experiment_real_gsa.avi');
time = toc

% For showing evolution of clusters where added and removed.
% save('./iskmeans_real_life1_gsa_show_evolution_of_clusters.mat', 'iCVs_ff', 'iCVs1_ff', 'ClusterIndices', 'labels', 'added', 'removed');


[ClusterIndices, labels] = ReconstructLabels(ClusterIndices, labels, false);
rand_index(ClusterIndices, labels,'adjusted')

FindNMI(ClusterIndices, labels)