% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h14_dq', 'var') && h14_dq > 0
    figure(h14_dq);
else
    h14_dq = figure;
end

ax(1) = subplot(2, 1, 1); hold on
plot(feature_full.t, dq_opt_plot_merge);
plot(t_recon_plot_merge, dq_recon_plot_merge, '.', 'LineWidth', 1);
title(['Trajectory and Reconstruction (RMSE (blended, windowed) : ' ...
    num2str(rmse_report.blended_rmse_belowThreshold_dq) ', ' num2str(rmse_report.windowed_rmse_belowThreshold_dq_mean) ' \pm ' num2str(rmse_report.windowed_rmse_belowThreshold_dq_std) ')']);
ylabel('Joint angle [rad]','Interpreter','Latex');
plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0)), ...
    zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0))) - 2, 'rx');
plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1)), ...
    zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1))) - 2, 'yx');
plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1)), ...
    zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1))) - 2, 'gx');

% switch run_mode
%     case 'win'
%         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2, -2);
% end

hold on
headers = {'Ankle $$q_{obs}$$', 'Knee $$q_{obs}$$', 'Hip $$q_{obs}$$', 'Ankle $$\hat{q}_{obs}$$', 'Knee $$\hat{q}_{obs}$$', 'Hip $$\hat{q}_{obs}$$', 'Pass = 0', 'Pass = 1', 'Pass \> 1'};
legend(headers,'Interpreter','Latex');

ax(2) = subplot(2, 1, 2); hold on
% title('Recovered Contribution Weights');
bar(feature_full.t, avgWeightArray_belowThres, 'stacked');
shading flat
ylim([0 1])
ylabel('Recovered $$\hat{c}$$','Interpreter','Latex');
xlabel('Time [s]');
legend(cost_function_names);

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);