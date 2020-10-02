% % % % plot data
% % % [~, uniqueInd] = unique(t_rmse);
% % % ax = [];
% % % 
% % % if exist('h1', 'var') && h1 > 0
% % %     figure(h1);
% % %     currXLim = xlim;
% % %     xlim([currXLim(1) feature_full.t(end)]);
% % % else
% % %     h1 = figure;
% % % end
% % % 
% % % ax(1) = subplot(2, 1, 1); hold on
% % % plot(feature_full.t, q_opt_plot_merge);
% % % plot(t_recon_plot_merge, q_recon_plot_merge, '.', 'LineWidth', 1);
% % % title(['RMSE Mean: ' num2str(rmse_report.meanRMSE) ', STD: ' num2str(rmse_report.stdRMSE), ', Range: ' num2str(rmse_report.rangeRMSE)]);
% % % 
% % % switch run_mode
% % %     case 'sim'
% % %         plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)));
% % %         
% % %     case 'win'
% % %         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
% % % end
% % % 
% % % ax(2) = subplot(2, 1, 2); hold on
% % % title('Contribution ratio (by sliding)');
% % % bar(t_rmse(uniqueInd), cAll_plot(uniqueInd, :), 'stacked'); 
% % % % bar(t_rmse(pts_underConstThres), cAll_plot(pts_underConstThres, :), 'stacked'); 
% % % % plot(t_rmse(pts_underConstThres), zeros(size(find(pts_underConstThres))), 'gx');
% % % 
% % % % if exist('q_jh_plot', 'var')
% % % %     threshold = 0.8;
% % % %     ptsPassThreshold = q_jh_plot >= threshold;
% % % %     ptsUnderThreshold = q_jh_plot < threshold;
% % % %     plot(t_rmse(ptsPassThreshold), zeros(size(find(ptsPassThreshold))) -0.05, 'gx');
% % % %     plot(t_rmse(ptsUnderThreshold), zeros(size(find(ptsUnderThreshold)))-0.05, 'rx');
% % % % end
% % % 
% % % shading flat
% % % legend(cost_function_names);
% % % 
% % % switch run_mode
% % %     case 'sim'
% % %         
% % %     case 'win'
% % %         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', [0 0], 1.05, -0.05);
% % % end
% % % % ylim([-10 110])
% % % 
% % % linkaxes(ax, 'x');
% % % % xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);

% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h1', 'var') && h1 > 0
    figure(h1);
else
    h1 = figure;
end

jointsToPlot = [4,5,6,9,10,11]; % back_jFB and right side arm and leg joints

ax(1) = subplot(2, 1, 1); hold on; grid on;
plot(feature_full.t, q_opt_plot_merge(jointsToPlot,:),'LineWidth',2);
legend(jointNames,'Location','northwest');
plot(t_recon_plot_merge, q_recon_plot_merge(:,jointsToPlot), '.', 'LineWidth', 1);
title(['Trajectory and Reconstruction (RMSE (blended, windowed) : ' ...
    num2str(rmse_report.blended_rmse_belowThreshold) ', ' num2str(rmse_report.windowed_rmse_belowThreshold_mean) ' \pm ' num2str(rmse_report.windowed_rmse_belowThreshold_std) ')']);
ylabel('Joint angle [rad]','Interpreter','Latex');
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0,  1.95, -1.95);

% switch run_mode
%     case 'win'
%         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2, -2);
% end

hold on
headers = {'Ankle $$q_{obs}$$', 'Knee $$q_{obs}$$', 'Hip $$q_{obs}$$', 'Ankle $$\hat{q}_{obs}$$', 'Knee $$\hat{q}_{obs}$$', 'Hip $$\hat{q}_{obs}$$', 'Pass = 0', 'Pass = 1', 'Pass \> 1'};
% legend(headers,'Interpreter','Latex');

ax(2) = subplot(2, 1, 2); hold on
% title('Recovered Contribution Weights');
barObj = bar(feature_full.t, avgWeightArray_belowThres(:, cost_function_names_sorted_ind), 'stacked');
shading flat
for i = 1:numel(cost_function_names_sorted_ind)
    barObj(i).FaceColor = colorVec(i,:);
end
ylim([0 1])
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 1.05, -0.05);
ylabel('Recovered $$\hat{c}$$','Interpreter','Latex');
xlabel('Time [s]');
legend(cost_function_names_sorted,'Location','northwest');

linkaxes(ax, 'x');
xlim([0 8]);
