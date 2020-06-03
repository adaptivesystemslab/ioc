% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

[sortVal, sortInd] = sort(cost_function_names);

allDof_plot_new = zeros(size(allDof_plot));
minDof_plot_new = zeros(size(minDof_plot));
for indla = 1:length(sortInd)
    allDof_plot_new(allDof_plot == indla) = sortInd(indla);
    minDof_plot_new(minDof_plot == indla) = sortInd(indla);
end

for indla = 1:length(cost_function_names)
    currCol = allDof_plot_new(:, indla);
    allDof_tallyWeights(indla) = length(find(currCol > 0));
    
    currCol = allDof_plot_new(:, indla);
    minDof_tallyWeights(indla) = length(find(minDof_plot > 0));
end


h2_2 = figure;
bar(tallyWeights); hold on
plot(t_rmse, minDof_plot_new, 'bo');
plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, max(max(allDof_plot)), 0);


% 
% ax(2) = subplot(2, 1, 2); hold on
% % title('Recovered Contribution Weights');
% bar(feature_full.t, avgWeightArray_belowThres(:, sortInd), 'stacked');

title('Selected pivot (b) over all the DOFs selected (r), alphabetrical order');

% figure; plot(minRmseInd_plot, 'bo', 'LineWidth', 2);
% figure; plot(minRmseIndArray, 'bo', 'LineWidth', 2);