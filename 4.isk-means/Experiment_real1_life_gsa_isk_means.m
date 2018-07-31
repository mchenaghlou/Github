clear all
filename = './dataset_real_gsa_isk_means.mat';
load(filename);
dim = size(gsa_isk_means,2) - 1;
labels = gsa_isk_means(:, end);
dataset = gsa_isk_means(:, 1:end-1);


% dataset = dataset(1:30000, :);
% labels = labels(1:30000, :);
% dataset(19551:24007,:) = [,];

% f=reduced_subset;
% inds = randi(size(f,1), 1,1500);
% ff =f(inds, 1:end);
% D = pdist(ff);
% D = squareform(D);
% clear ff
% vativat(D, 'meteorite landing');
tic
% [1:15000,16007:17007,24008:25008,34009:40000]
[ V, ClusterIndices, iCVs_ff, iCVs1_ff, added, removed] = isKmeans(dataset, 0, 0.99, 0.99, 0, ...
    './4.isk-means/Experiment_real_gsa.avi');
time = toc
% For showing evolution of clusters where added and removed.
save('./4.isk-means/PaperExperiments/iskmeans_real_life1_gsa_show_evolution_of_clusters.mat', 'iCVs_ff', 'iCVs1_ff', 'ClusterIndices', 'labels', 'added', 'removed');
% 
% figure('rend','painters','pos',[300 600 400 300])
% hold on
% h1 = plot(iCVs_ff)
% h2 = plot(iCVs1_ff)
% 
% legend('iXB_{\lambda} isk-means','iXB_{\lambda} cisk-means')
% 
% set(h1,'linewidth',2);
% set(h2,'linewidth',2);
% 
% xlim([0 63014])
% ylim([0 20])
% 
% % xt = get(gca, 'XTick');
% set(gca, 'FontSize', 14)
% set(gca, 'FontWeight', 'bold')
% 
% % yt = get(gca, 'YTick');
% set(gca, 'FontSize', 14)
% set(gca, 'FontWeight', 'bold')
% 
% xlabel('Time', 'FontSize',10,'FontWeight','bold')
% ylabel('iXB_{\lambda} values', 'FontSize',10,'FontWeight','bold')
% 
% unique(ClusterIndices)
% % labels = reduced_subset([1:15000,34009:40000], end );
% unique(labels)
% 
% XB_batch(dataset, V)


[ClusterIndices, labels] = ReconstructLabels(ClusterIndices, labels, false);
rand_index(ClusterIndices, labels,'adjusted')

FindNMI(ClusterIndices, labels)