% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

h12 = figure;
ax(1) = subplot(3, 1, 1); hold on
plot(feature_full.t, q_opt_plot_merge);
plot(t_recon_plot_merge, q_recon_plot_merge, '.', 'LineWidth', 1);
title(['RMSE Mean: ' num2str(rmse_report.meanRMSE) ', STD: ' num2str(rmse_report.stdRMSE), ', Range: ' num2str(rmse_report.rangeRMSE)]);

switch run_mode
    case 'sim'
        plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)));
        
    case 'win'
        plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
end

ax(2) = subplot(3, 1, 2); hold on
title('Recovered contribution weight');
bar(feature_full.t, avgWeightArray_belowThres, 'stacked'); 
shading flat
ylim([0 1])

ax(3) = subplot(3, 1, 3); hold on
plotType = 'ratio';
plotCalc_metricCalc;
title('Recovered contribution ratio');
bar(feature_full.t, avgWeightArray_belowThres, 'stacked');
shading flat
ylim([0 1])
legend(cost_function_names);

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);