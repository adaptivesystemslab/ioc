% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h1_2', 'var') && h1_2 > 0
    figure(h1_2);
    currXLim = xlim;
    xlim([currXLim(1) feature_full.t(end)]);
else
    h1_2 = figure;
end

ax(1) = subplot(3, 1, 1); hold on
plot(feature_full.t, q_opt_plot_merge);
plot(t_recon_plot_merge, q_recon_plot_merge, '.', 'LineWidth', 1);
title(['(sliding) RMSE Mean: ' num2str(rmse_report.meanRMSE) ', STD: ' num2str(rmse_report.stdRMSE), ', Range: ' num2str(rmse_report.rangeRMSE)]);

switch run_mode
    case 'sim'
        plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)));
        
    case 'win'
        plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
end

ax(2) = subplot(3, 1, 2); hold on
title('Contribution ratio (by sliding)');
bar(t_rmse(uniqueInd), cAll_plot(uniqueInd, :), 'stacked'); 
% bar(t_rmse(pts_underConstThres), cAll_plot(pts_underConstThres, :), 'stacked'); 
% plot(t_rmse(pts_underConstThres), zeros(size(find(pts_underConstThres))), 'gx');

% if exist('q_jh_plot', 'var')
%     threshold = 0.8;
%     ptsPassThreshold = q_jh_plot >= threshold;
%     ptsUnderThreshold = q_jh_plot < threshold;
%     plot(t_rmse(ptsPassThreshold), zeros(size(find(ptsPassThreshold))) -0.05, 'gx');
%     plot(t_rmse(ptsUnderThreshold), zeros(size(find(ptsUnderThreshold)))-0.05, 'rx');
% end

shading flat
legend(cost_function_names);

switch run_mode
    case 'sim'
        
    case 'win'
        plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', [0 0], 1.05, -0.05);
end
% ylim([-10 110])



ax(3) = subplot(3, 1, 3); hold on
title('Resnorm (by sliding)');
bar(t_rmse(uniqueInd), resnorm_min(uniqueInd, :), 'stacked'); 
hold on
ylimSample = ylim;
plot(xlim, ones(2, 1)*constThres, 'k');
ylim(ylimSample);

shading flat


linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);

