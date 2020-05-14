% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

h7 = figure;
ax(1) = subplot(2, 1, 1); hold on
title('Recovered contribution ratio');
bar(t_rmse(uniqueInd), ratioRMSE_plot(uniqueInd, :), 'stacked'); 

plot(t_rmse(pts_deny), zeros(size(find(pts_deny)))-0.05, 'rx');
plot(t_rmse(pts_allow), zeros(size(find(pts_allow)))  -0.05, 'gx');

switch run_mode
    case 'sim'
%         plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)));
        
    case 'win'
        plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 105, -5);
end

ylim([-10 110])



shading flat
% 
% switch run_mode
%     case 'sim'
%         plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)));
%         
%     case 'win'
%         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
% end

ax(2) = subplot(2, 1, 2); hold on
area(t_rmse(uniqueInd), resnormAll_lsqlin_const_RMSE_plot(uniqueInd, :), 'FaceColor', [0 0 1]); 
% area(t_rmse(uniqueInd), resnormAll_lsqlin_unconst_RMSE_plot(uniqueInd, :), 'FaceColor', [0 1 0]); 
area(feature_full.t, resnromAll_lsqlin_const_minRMSE_array(:, :), 'FaceColor', [0 1 0]); 

title(['Const (b) resnorm (mu=' num2str(meanconstResnrom) '), const avg (g) resnorm (mu=' num2str(meanunconstResnorm) ') and threshold (k)']);
plot(xlim, [constThres constThres], 'k');

% plot(t_rmse(pts_allow), zeros(size(find(pts_allow)))  -0.05, 'gx');
% plot(t_rmse(pts_deny), zeros(size(find(pts_deny)))    -0.05, 'rx');

% plot(t_rmse(pts_aboveConstThres), zeros(size(find(pts_aboveConstThres)))  -0.0001, 'rx');
plot(t_rmse(pts_underConstThres), zeros(size(find(pts_underConstThres)))  -0.0001, 'mx');

% plot(t_rmse(pts_aboveDiffThres), zeros(size(find(pts_aboveDiffThres)))  -0.0003, 'rx');
plot(feature_full.t(pts_underConstThres_averageWindow), zeros(size(find(pts_underConstThres_averageWindow))) + 0, 'rx');

ylim([-0.0005 constThres*5]);


% legend('show');

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);