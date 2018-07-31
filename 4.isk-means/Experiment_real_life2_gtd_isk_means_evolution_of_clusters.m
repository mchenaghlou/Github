clear all
load('./gtd_for_iskmeans_europe.mat');
load('./iskmeans_real_life2_gtd_show_evolution_of_clusters.mat');


%% 1. Plot iXB_lambda values
f = figure('pos',[1000 200 900 400]);
% p = uipanel('Parent',f,'BorderType','none');
% subplot(2,1,1,'Parent',p)
deviate = 0.2;
y_limit = 4;
plot_icvs_ffs( iCVs_ff,iCVs1_ff, added, removed, deviate, y_limit)
% xlim([0 5324])
%%%%%%%%
begin = 4500;
ending = 4800;
% figure
% plot(dataset(4500:4800,1), dataset(4500:4800, 2), '.')
rect_height = 1.2;
rectangle('Position',[begin-10 0 ending-begin+10 rect_height])
line([begin-10, 1175],[rect_height 2.62], 'color', 'k','HandleVisibility','off')
line([ending, 3230],[rect_height 2.62], 'color', 'k','HandleVisibility','off')
axes('position',[.3 .65 .3 .25])
box on % put box around new pair of axes
hold on
h1 = plot(begin:ending, iCVs_ff(begin:ending))
h2 = plot(begin:ending, iCVs1_ff((begin:ending)))
set(h1,'linewidth',1.5);
set(h2,'linewidth',1.5);
set(gca, 'FontSize', 12)
set(gca, 'FontWeight', 'bold')
period_added = added(added <= ending & added >= begin)
period_removed = removed(removed <= ending & removed >= begin)
% for i = 1:length(period_added)
%     temp1 = period_added(i);
%     plot(temp1, iCVs1_ff(temp1)+ 0.05, '+r')    
%     text(temp1, iCVs1_ff(temp1)+ 0.1,num2str(ClusterIndices(period_added(i)+1)));
% end
plot(period_added, iCVs1_ff(period_added)+ 0.05, '+r')
plot(period_removed, iCVs_ff(period_removed)+ 0.05, 'vk')
axis tight
xlim([begin ending]);
xticks(begin:100:ending)
ylim([0.05 0.85]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. Plot Cluster Evolution of Terror Attacks
figure
hold on

for i = 1:length(unique(ClusterIndices))
%     if sum(ClusterIndices == i) > 0
        curr_i = unique(ClusterIndices);
        curr_i = curr_i(i);
        x_values = find(ClusterIndices' == curr_i)
        y_values = curr_i .* ones(1,length(x_values));
        plot(x_values, y_values , '.')
%     end
end
xlim([0 5324])
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

% yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

ylabel('Cluster Index', 'FontSize',14,'FontWeight','bold')
xlabel('Time', 'FontSize',14,'FontWeight','bold')

%% 3. Plot Clustering
figure
hold on
for i = 1:length(unique(ClusterIndices))
    plot(dataset(ClusterIndices == i, 1), dataset(ClusterIndices == i, 2), '.')
end
for i = 1:size(V, 1)
    plot(V(i, 1), V(i, 2), '*k', 'markersize', 12,'linewidth', 2)
    text(V(i, 1), V(i, 2) + 1.5,num2str(ClusterLabelDic(ClusterLabelDic(:, 1) == i, 2)), 'FontWeight', 'bold');
end
 


set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

% yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

ylabel('Longitude', 'FontSize',14,'FontWeight','bold')
xlabel('Latitude', 'FontSize',14,'FontWeight','bold')

%% 4. Plot Cluster Evolution

% for j = 1:length(cluster_evolution_snapshots.prototypes)
%     figure
%     hold on
% 
%     V_temp = cluster_evolution_snapshots.prototypes{j};
%     datas = cluster_evolution_snapshots.data{j};
%     for k = 1:length(unique(ClusterIndices(1:size(datas, 1))))
%         plot(datas(ClusterIndices(1:size(datas, 1)) == k, 1), datas(ClusterIndices(1:size(datas, 1)) == k, 2), '.')
%     end
%     for i = 1:size(V_temp, 1)
%         plot(V_temp(i, 1), V_temp(i, 2), '*k', 'markersize', 12,'linewidth', 2)
%         text(V_temp(i, 1), V_temp(i, 2) + 2,num2str(i), 'FontWeight', 'bold');
%     end
% end





