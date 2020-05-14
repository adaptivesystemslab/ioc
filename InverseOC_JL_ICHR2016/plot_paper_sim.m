% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h13', 'var') && h13 > 0
    figure(h13);
else
    h13 = figure;
end
 
% ax(1) = subplot(2, 1, 1); hold on
% plot(feature_full.t, q_opt_plot_merge);
% % plot(t_recon_plot_merge, q_recon_plot_merge, '.', 'LineWidth', 1);
% % plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0)), ...
% %     zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0))) - 2, 'rx');
% % plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1)), ...
% %     zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1))) - 2, 'yx');
% % plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1)), ...
% %     zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1))) - 2, 'gx');
% switch run_mode
%     case 'win'
%         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2, -2);
% end
% % title(['Trajectory and Reconstruction (RMSE: ' num2str(rmse_report.meanRMSE) ' \pm ' num2str(rmse_report.stdRMSE) ')']);
% ylabel('Joint angle [rad]','Interpreter','Latex');
% % headers = {'Ankle $$q_{obs}$$', 'Knee $$q_{obs}$$', 'Hip $$q_{obs}$$', 'Ankle $$\hat{q}_{obs}$$', 'Knee $$\hat{q}_{obs}$$', 'Hip $$\hat{q}_{obs}$$'};
% headers = {'Ankle $$q_{obs}$$', 'Knee $$q_{obs}$$', 'Hip $$q_{obs}$$'};
% legend(headers,'Interpreter','Latex');



ax(1) = subplot(2, 1,1); hold on
bar(feature_full.t', ccost_full_plot', 'stacked'); 
shading flat
ylim([0 1])
legend(cost_function_names);
% title('Simulated Contribution Weights');
ylabel('Generated $${c}$$','Interpreter','Latex');
% title('Recovered Contribution Weights');

ax(2) = subplot(2, 1, 2); hold on
bar(feature_full.t, avgWeightArray_belowThres, 'stacked');
shading flat
ylim([0 1])
ylabel('Recovered $$\hat{c}$$','Interpreter','Latex');
xlabel('Time [s]');

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);