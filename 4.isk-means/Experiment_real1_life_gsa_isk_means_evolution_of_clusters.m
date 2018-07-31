clear all
load('./iskmeans_real_life1_gsa_show_evolution_of_clusters.mat');
labels = ChangeLabelsFrom1ToN( labels );


f = figure('pos',[1000 200 900 400]);

deviate = 0.5;
y_limit = 25;
plot_icvs_ffs( iCVs_ff,iCVs1_ff, added, removed, deviate, y_limit)
begin = 16490;
ending = 24800;

rectangle('Position',[begin-10 0 ending-begin 5])
line([begin-10, 15500],[5 12.75], 'color', 'k','HandleVisibility','off')
line([ending, 35750],[5 12.75], 'color', 'k','HandleVisibility','off')

axes('position',[.32 .55 .25 .3])

box on % put box around new pair of axes
hold on

h1 = plot(begin:ending, iCVs_ff(begin:ending))
h2 = plot(begin:ending, iCVs1_ff((begin:ending)))


set(h1,'linewidth',2);
set(h2,'linewidth',2);
set(gca, 'FontSize', 12)
set(gca, 'FontWeight', 'bold')

period_added = added(added <= ending & added >= begin)
period_removed = removed(removed <= ending & removed >= begin)
plot(period_added, iCVs1_ff(period_added)+ deviate, '+r')
plot(period_removed, iCVs_ff(period_removed)+ deviate, 'vk')

axis tight
xlim([begin ending]);
xticks([begin+10:2000:ending])
ylim([0 5]);
% subplot(2,1,2,'Parent',p)
% plot_icvs_ffs( iCVs_ff,iCVs1_ff, added, removed, 595, 640)
