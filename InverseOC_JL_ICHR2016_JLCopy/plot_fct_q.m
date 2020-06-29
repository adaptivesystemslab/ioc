function h14 = plot_fct_q(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, dq_obs, dq_doc, ddq_obs, ddq_doc, ...
    rmsEntryToUse, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd)

ax = [];
h14 = figure;

rmsUse = 1:size(t_doc, 2);
% rmsUse = rmsUse;

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


ax(2) = subplot(3, 1, 2); hold on; grid on;
for i = 1:numel(jointsToPlot)
    plot(t_obs, dq_obs(jointsToPlot(i),:),'LineWidth',2,'Color',colorVec2(i,:));
    % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
    plot(t_doc(rmsUse), dq_doc(jointsToPlot(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
end
legend(jointNames,'Location','northwest');
title(['dq']);
ylabel('Joint velo [rad/s]','Interpreter','Latex');
% plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);


ax(3) = subplot(3, 1, 3); hold on;  grid on;
for i = 1:numel(jointsToPlot)
    plot(t_obs, ddq_obs(jointsToPlot(i),:),'LineWidth',2,'Color',colorVec2(i,:));
    % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
    plot(t_doc(rmsUse), ddq_doc(jointsToPlot(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
end
legend(jointNames,'Location','northwest');
title(['ddq']);
ylabel('Joint accel [rad/$$s^2$$]','Interpreter','Latex');
% plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);

linkaxes(ax, 'x');
xlim([t_obs(1) t_obs(end)]);