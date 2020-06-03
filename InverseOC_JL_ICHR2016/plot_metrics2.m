% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

h2_2 = figure;
ax(2) = subplot(2, 2, 1); hold on
plot(t_rmse, rmseMin_plot, 'k', 'LineWidth', 2); ylimGet = ylim; % plot first time to generate the ylim data
% barData = accumarray(t_rmse,rmseAll_plot,[],@(x) {hist(x,ageValues)});
% bar(t_rmse, rmseAll_plot, 'BarWidth', 2, 'BarLayout', 'grouped'); ylim(ylimGet); shading flat
plot(t_rmse, rmseAll_plot); ylim(ylimGet);
plot(t_rmse, rmseMin_plot, 'k', 'LineWidth', 2);
title('rmse, norm to length & dof');

ax(3) = subplot(2, 2, 2); hold on
t_rmseStack = repmat(t_rmse, 1, size(dofsUsedInIOC_plot, 2));
dofsUsedInIOC_plotStack = reshape(dofsUsedInIOC_plot, 1, size(dofsUsedInIOC_plot, 1)*size(dofsUsedInIOC_plot, 2));
plot(t_rmseStack(dofsUsedInIOC_plotStack > 0), dofsUsedInIOC_plotStack(dofsUsedInIOC_plotStack > 0), 'kx', 'LineWidth', 2); % plot first time to generate the ylim data
% barData = accumarray(t_rmse,rmseAll_plot,[],@(x) {hist(x,ageValues)});
% bar(t_rmse, rmseAll_plot, 'BarWidth', 2, 'BarLayout', 'grouped'); ylim(ylimGet); shading flat
plot(t_rmse, minRmseInd_plot, 'bo', 'LineWidth', 2);
title('Pivots tested (k), with selected pivot (b)');
ylim([0.5 length(ccost_win)+0.5]);

ax(5) = subplot(2, 2, 3); hold on
plot(t_rmse, resnormMin_plot, 'k', 'LineWidth', 2); ylimGet = ylim;
plot(t_rmse, resnormAll_plot); ylim(ylimGet);
plot(t_rmse, resnormMin_plot, 'k', 'LineWidth', 2);
plot(xlim, ones(2, 1)*constThres, 'k');
title('resnorm');

% ax(6) = subplot(2, 2, 4); hold on
% plot(t_rmse, condMin_plot, 'k', 'LineWidth', 2); ylimGet = ylim;
% plot(t_rmse, condAll_plot); ylim(ylimGet);
% plot(t_rmse, condMin_plot, 'k', 'LineWidth', 2);
% title('cond');

ax(5) = subplot(2, 2, 4); hold on
plot(t_rmse, resnormAll_lsqlin_const_plot, 'k', 'LineWidth', 2); 
plot(xlim, ones(2, 1)*constThres, 'k');
title('resnorm lsqlin');
ylim([0 constThres*10]);

linkaxes(ax, 'x');
xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);