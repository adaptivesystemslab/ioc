% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

[sortVal, sortInd] = sort(cost_function_names);

allDof_plot_new = zeros(size(allDof_plot));
minDof_plot_new = zeros(size(minDof_plot));
for indla = 1:length(sortInd)
    allDof_plot_new(allDof_plot == sortInd(indla)) = indla;
    minDof_plot_new(minDof_plot == sortInd(indla)) = indla;
    
    allDof_tallyWeights(indla) = length(find(allDof_plot(:, indla) > 0));
    minDof_tallyWeights(indla) = length(find(minDof_plot == indla));
end

% h2_1 = figure; % orig
% plot(t_rmse, allDof_plot, 'rx'); hold on
% plot(t_rmse, minDof_plot, 'bo');
% ylim([0.5 length(cost_function_names)+0.5]);
% set(gca,'YTickLabel',cost_function_names);

h2_2 = figure;
plot(t_rmse, allDof_plot_new, 'bx'); hold on
plot(t_rmse, minDof_plot_new, 'ro');
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, max(max(allDof_plot)), 0);
ylim([0.5 length(cost_function_names)+0.5]);
getYLab = get(gca, 'YTickLabel');
ylab = str2num(getYLab);

labUpdater = {}; labUpdaterCounter = 0;
for ind_ylab = 1:length(ylab)
    if mod(ylab(ind_ylab), 1) == 0
        labUpdaterCounter = labUpdaterCounter + 1;
        labUpdater{ind_ylab} = sortVal{labUpdaterCounter};
    else
        labUpdater{ind_ylab} = '';
    end
end

set(gca,'YTickLabel',labUpdater);
title('Selected pivot (r) over all the DOFs selected (b), alphabetrical order');

% h2_3 = figure;
% bar(allDof_tallyWeights); hold on
% set(gca,'XTickLabel',cost_function_names);
% hh = rotateXLabels(gca, 90);

dofcomb = [allDof_tallyWeights; minDof_tallyWeights];
h2_4 = figure;
bar(dofcomb(:, sortInd)'); hold on
set(gca,'XTickLabel',sortVal);
hh = rotateXLabels(gca, 90);
title('Selected pivot (r) over all the DOFs selected (b), alphabetrical order');

% 
% ax(2) = subplot(2, 1, 2); hold on
% % title('Recovered Contribution Weights');
% bar(feature_full.t, avgWeightArray_belowThres(:, sortInd), 'stacked');



% figure; plot(minRmseInd_plot, 'bo', 'LineWidth', 2);
% figure; plot(minRmseIndArray, 'bo', 'LineWidth', 2);