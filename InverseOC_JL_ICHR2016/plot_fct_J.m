function h14 = plot_cost_fct(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, avgRatioArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
    rmsEntryToUse, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd)

ax = [];
h14 = figure;

rmsUse = 1:size(t_doc, 2);
% rmsUse = find(rmsEntryToUse);

ax(1) = subplot(3, 1, 1); hold on; grid on;
for i = 1:numel(jointsToPlot)
    plot(t_obs, q_obs(jointsToPlot(i),:),'LineWidth',2,'Color',colorVec2(i,:));
    % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
    plot(t_doc(rmsUse), q_doc(jointsToPlot(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
end
legend(jointNames,'Location','northwest');
title(['Trajectory and Reconstruction (RMSE: ' num2str(rmseMean) ' \pm ' num2str(rmseStd) ')']);
ylabel('Joint angle [rad]','Interpreter','Latex');
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);
ylim([-2.05, 2.05]);


ax(2) = subplot(3, 1, 2); hold on;
% title('Recovered Contribution Weights');
barObj = bar(t_obs, avgRatioArray_belowThres(:, cost_function_names_sorted_ind), 'stacked');
shading flat
for i = 1:numel(cost_function_names_sorted_ind)
    barObj(i).FaceColor = colorVec(i,:);
end
ylim([0 1])
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 1.05, -0.05);
ylabel('Recovered $$\hat{J}$$ ratio','Interpreter','Latex');
% xlabel('Time [s]');
legend(cost_function_names_sorted,'Location','northwest');


ax(3) = subplot(3, 1, 3); hold on; grid on;
area(t_obs, resnromAll_lsqlin_const_minRMSE_array_belowThres(:, :), 'FaceColor', [0 0 1]); 
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
la = ylim;
hold on
xlabel('Time [s]');
ylabel('Averaged Residual Norm');
title([' RESNORM  : ' num2str(resnormMean) ' \pm ' num2str(resnormStd)]);
ylim([0 la(2)]);

linkaxes(ax, 'x');
xlim([t_obs(1) t_obs(end)]);