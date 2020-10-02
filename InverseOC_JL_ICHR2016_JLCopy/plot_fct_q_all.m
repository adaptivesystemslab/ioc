function h14 = plot_fct_q_all(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, dq_obs, dq_doc, ddq_obs, ddq_doc, ...
    rmsEntryToUse, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd)

ax = [];
h14 = figure;

rmsUse = 1:size(t_doc, 2);
% rmsUse = rmsUse;

% torsoJointsPr = 1:2;  torsoJointNamesPr = {'pX','pZ'};
% torsoJointsRot = 3:4;  torsoJointNamesRot = {'rY','back-FB'};
armJointsR = 1:3;    armJointNamesR = {'sho-x', 'sho-y', 'sho-z'};
armJointsL = 4;    armJointNamesL = {'elb-x'};
legJointsR = 5:7;   legJointNamesR = {'wri-x', 'wri-y', 'wri-z'};
% legJointsL = 12:14;   legJointNamesL = {'hip-L','knee-L','ankle-L'};

ax(1) = subplot(3, 1, 1); hold on; grid on;
for i = 1:numel(armJointsR)
    plot(t_obs, q_obs(armJointsR(i),:),'LineWidth',2,'Color',colorVec2(i,:));
    % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
    plot(t_doc(rmsUse), q_doc(armJointsR(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
end
legend(armJointNamesR,'Location','northwest');
title(['shoulder']);
ylabel('Joint angle [rad]','Interpreter','Latex');
% plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);
% ylim([-2.05, 2.05]);


ax(2) = subplot(3, 1, 2); hold on; grid on;
for i = 1:numel(armJointsL)
    plot(t_obs, q_obs(armJointsL(i),:),'LineWidth',2,'Color',colorVec2(i,:));
    % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
    plot(t_doc(rmsUse), q_doc(armJointsL(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
end
legend(armJointNamesL,'Location','northwest');
title(['elbow']);
ylabel('Joint angle [rad]','Interpreter','Latex');
% plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);
% ylim([-2.05, 2.05]);


ax(3) = subplot(3, 1, 3); hold on;  grid on;
for i = 1:numel(legJointsR)
    plot(t_obs, q_obs(legJointsR(i),:),'LineWidth',2,'Color',colorVec2(i,:));
    % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
    plot(t_doc(rmsUse), q_doc(legJointsR(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
end
legend(legJointNamesR,'Location','northwest');
title(['wrist']);
ylabel('Joint angle [rad]','Interpreter','Latex');
% plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);
% ylim([-2.05, 2.05]);


% ax(6) = subplot(3, 2, 6); hold on;  grid on;
% for i = 1:numel(legJointsL)
%     plot(t_obs, q_obs(legJointsL(i),:),'LineWidth',2,'Color',colorVec2(i,:));
%     % plot(t_recon_plot_merge, q_recon_plot_merge(jointsToPlot,:), '.', 'LineWidth', 1);
%     plot(t_doc(rmsUse), q_doc(legJointsL(i), rmsUse), '.', 'LineWidth', 1,'Color',colorVec2(i,:),'HandleVisibility','off');
% end
% legend(legJointNamesL,'Location','northwest');
% title(['leg joints left']);
% ylabel('Joint angle [rad]','Interpreter','Latex');
% plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2.05, -2.05);
% ylim([-2.05, 2.05]);


linkaxes(ax, 'xy');
xlim([t_obs(1) t_obs(end)]);