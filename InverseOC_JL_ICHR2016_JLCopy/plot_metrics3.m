% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

h2_3 = figure;
ax(1) = subplot(211);
dofsUsedInIOC_plotStack = reshape(dofsUsedInIOC_plot, 1, size(dofsUsedInIOC_plot, 1)*size(dofsUsedInIOC_plot, 2));
hist(dofsUsedInIOC_plotStack(dofsUsedInIOC_plot > 0), 1:length(cost_function_names)); % plot first time to generate the ylim data
% barData = accumarray(t_rmse,rmseAll_plot,[],@(x) {hist(x,ageValues)});
% bar(t_rmse, rmseAll_plot, 'BarWidth', 2, 'BarLayout', 'grouped'); ylim(ylimGet); shading flat
title('Number of times the pivots was selected for testing');

ax(2) = subplot(212);
hist(minRmseInd_plot, 1:length(cost_function_names));
title('Number of times the pivots was selected as minRMSE');

linkaxes(ax, 'x');
xlim([0.5 length(cost_function_names) + 0.5]);

subplot(211);
set(gca,'XTickLabel',cost_function_names);
hh = rotateXLabels( gca(), 45 );

subplot(212);
set(gca,'XTickLabel',cost_function_names);
hh = rotateXLabels( gca(), 45 );