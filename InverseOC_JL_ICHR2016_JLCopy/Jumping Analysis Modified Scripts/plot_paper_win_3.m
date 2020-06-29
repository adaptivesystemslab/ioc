% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h14', 'var') && h14 > 0
    figure(h14);
else
    h14 = figure;
end

jointNames = {'back-FB','shldr','elbow','hip','knee','ankle'};
jointsToPlot = [4,5,6,9,10,11]; % back_jFB and right side arm and leg joints

ax(1) = subplot(3, 1, 1); hold on; grid on;
plot(feature_full.t, q_opt_plot_merge(jointsToPlot,:),'LineWidth',2);
legend(jointNames,'Location','northwest');
% plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
plot(feature_full.t(find(rmsEntryToUse)), q_blended(jointsToPlot, find(rmsEntryToUse)), '--', 'LineWidth', 1);
title(['Trajectory and Reconstruction discarding belowthreshold values (RMSE (blended, windowed) : ' ...
    num2str(rmse_report.blended_rmse_belowThreshold) ', ' num2str(rmse_report.windowed_rmse_belowThreshold_mean) ' \pm ' num2str(rmse_report.windowed_rmse_belowThreshold_std) ')']);
ylabel('Joint angle [rad]','Interpreter','Latex');
plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0)), ...
    zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0))) - 2, 'rx');
plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1)), ...
    zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1))) - 2, 'yx');
plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1)), ...
    zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1))) - 2, 'gx');
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);

% switch run_mode
%     case 'win'
%         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2, -2);
% end

hold on
headers = {'Ankle $$q_{obs}$$', 'Knee $$q_{obs}$$', 'Hip $$q_{obs}$$', 'Ankle $$\hat{q}_{obs}$$', 'Knee $$\hat{q}_{obs}$$', 'Hip $$\hat{q}_{obs}$$', 'Pass = 0', 'Pass = 1', 'Pass \> 1'};
% legend(headers,'Interpreter','Latex');



ax(2) = subplot(3, 1, 2); hold on;
% title('Recovered Contribution Weights');
barObj = bar(feature_full.t, avgWeightArray_belowThres(:, cost_function_names_sorted_ind), 'stacked');
shading flat
for i = 1:numel(cost_function_names_sorted_ind)
    barObj(i).FaceColor = colorVec(i,:);
end
ylim([0 1])
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 1.05, -0.05);
ylabel('Recovered $$\hat{c}$$','Interpreter','Latex');
% xlabel('Time [s]');
legend(cost_function_names_sorted,'Location','northwest');

ax(3) = subplot(3, 1, 3); hold on; grid on;
area(feature_full.t, resnromAll_lsqlin_const_minRMSE_array_belowThres(:, :), 'FaceColor', [0 0 1]); 
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
la = ylim;
hold on
% title('Averaged Residual Norm');
plot(xlim, [constThres constThres], 'k');
xlabel('Time [s]');
ylabel('Averaged Residual Norm');
title([' RESNORM  : ' ...
    num2str(rmse_report.resnorm_mean) ' \pm ' num2str(rmse_report.resnorm_std)]);
ylim([0 la(2)]);

linkaxes(ax, 'x');
xlim([0 8]);
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);