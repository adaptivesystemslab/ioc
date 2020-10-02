% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

h2_1 = figure;
ax(1) = subplot(2, 1, 1); hold on
bar(traj_load_lengthAdj.t, J_full, 'stacked'); hold on
title('Simulated contribution ratio');
shading flat
legend(cost_function_names);

ax(4) = subplot(2, 1, 2); hold on
title('Actual contribution ratio');
bar(t_rmse(uniqueInd), ratioRMSE_plot(uniqueInd, :), 'stacked');

shading flat
legend(cost_function_names);
ylim([-10 110])

linkaxes(ax, 'x');
xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);