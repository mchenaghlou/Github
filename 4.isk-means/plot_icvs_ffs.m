function [ output_args ] = plot_icvs_ffs(iCVs_ff, iCVs1_ff, added, removed, deviation, y_limit, begin, ending)
%PLOT_ICVS_FFS Summary of this function goes here
%   Detailed explanation goes here


% f = figure('rend','painters','pos',[300 300 800 250])
hold on
switch nargin
    case 6
        begin = 1;
        ending = length(iCVs_ff);
end
h1 = plot(begin:ending, iCVs_ff(begin:ending))
h2 = plot(begin:ending, iCVs1_ff((begin:ending)))




period_added = added(added <= ending & added >= begin)
period_removed = removed(removed <= ending & removed >= begin)
h3 = plot(period_added, iCVs1_ff(period_added)+ deviation, '+r')
h4 = plot(period_removed, iCVs_ff(period_removed)+ deviation, 'vk')

% legend('iXB_{\lambda} isk-means','iXB_{\lambda} cisk-means', 'Prototype Added', 'Prototype Removed', 'Location', 'northeast')
leg = {'iXB_{\lambda} Current-sKM','iXB_{\lambda} Control-sKM', 'Prototype Added', 'Prototype Removed'}
% columnlegend(2, leg, 'boxon') 
% legend(leg)
legend([h1 h2 h3 h4],leg);


set(h1,'linewidth',2);
set(h2,'linewidth',2);

xlim([begin ending])
ylim([0 y_limit])

% xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

% yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')

xlabel('Time', 'FontSize',14,'FontWeight','bold')
ylabel('iXB_{\lambda} values', 'FontSize',14,'FontWeight','bold')

end

