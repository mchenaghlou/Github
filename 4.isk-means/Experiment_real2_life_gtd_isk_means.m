clear all

%%%%% script to make europe dataset.
% load("./../Data/DataSets_realworld/gtd/gtd.mat");
% % dataset = gtd;
% gtd = unique(gtd,'rows','stable');
% % plot(dataset(:, 1), dataset(:, 2), '.')
% init_xs = -7.5 <= gtd(:, 1) & gtd(:, 1) <= -5.5;
% init_ys =  54 <= gtd(:, 2) & gtd(:, 2) <= 55;
% init_pop = gtd(init_xs & init_ys, :) 
% % plot(init_pop(:, 1), init_pop(:, 2), '.')
% xs =  35 <= gtd(:, 2) & gtd(:, 2) <= 70;
% ys = -10 <= gtd(:, 1) & gtd(:, 1) <= 35;
% xy = xs & ys;
% dataset = [init_pop(1:80, :); gtd(xy, :)];
% % plot(dataset(:, 1), dataset(:, 2), '.')
% save('./../Data/DataSets_realworld/gtd/gtd_for_iskmeans_europe.mat', 'dataset')
%%%%% script to make europe dataset.


load('./GTD_dataset.mat');


 [ V, ClusterIndices, iCVs_ff, iCVs1_ff, added, removed, ...
     cluster_evolution_snapshots, ClusterLabelDic] = isKmeans(dataset, 2, ...
     0.99, 0.99, 0, './Experiment_real_gtd.avi');
% load('./4.isk-means/PaperExperiments/iskmeans_real_life2_gtd_show_evolution_of_clusters.mat');

% For showing evolution of clusters where added and removed.
save('./iskmeans_real_life2_gtd_show_evolution_of_clusters.mat', 'V','iCVs_ff', 'iCVs1_ff', 'ClusterIndices', 'added', 'removed', 'cluster_evolution_snapshots', 'ClusterLabelDic');

