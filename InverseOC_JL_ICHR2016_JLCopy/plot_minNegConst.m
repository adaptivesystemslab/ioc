% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

h9 = figure;
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
area(t_rmse(uniqueInd), minWeightsUnconst(uniqueInd, :), 'FaceColor', [0 0 1]); 
title(['Min unconstrained value and threshold (k)']);
plot(xlim, [negConstThreshold negConstThreshold], 'k');

% plot(t_rmse(pts_allow), zeros(size(find(pts_allow)))  -0.05, 'gx');
% plot(t_rmse(pts_deny), zeros(size(find(pts_deny)))    -0.05, 'rx');

plot(t_rmse(lowestValInUnconst), zeros(size(find(lowestValInUnconst)))  -0.0001, 'gx');


% ylim([-0.0005 constThres*5]);


% legend('show');

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);