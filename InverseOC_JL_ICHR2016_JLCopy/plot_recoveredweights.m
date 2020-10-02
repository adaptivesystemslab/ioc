% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

h6 = figure;
ax(1) = subplot(2, 1, 1); hold on
title('Contribution ratio');
bar(t_rmse(uniqueInd), c_ioc_constrain_plot(uniqueInd, :)); 
shading flat
ylim([-0.2 1.5])

ax(5) = subplot(2, 1, 2); hold on
title('Contribution ratio');
bar(t_rmse(uniqueInd), c_ioc_unconstrain_plot(uniqueInd, :)); 

shading flat
ylim([-0.2 1.5])

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);