% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

upperBound = length(cost_function_names);
if length(cost_function_names) > 8
    upperBound = 8;
end

h3_2 = figure;
for ii = 5:upperBound
    ax(ii-4) = subplot(2, 2, ii-4); hold on
    bar(t_rmse(uniqueInd), ratioRMSE_plot(uniqueInd, ii));
    shading flat
    title(['Contrib ratio: ' cost_function_names{ii}]);
    
    switch run_mode
        case 'sim'
            plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)), 'k', 0, 100.5, -0.5);

        case 'win'
            plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 100.5, -0.5);
    end
end

if ~isempty(ax)
    linkaxes(ax, 'xy');
end
xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);
ylim([-10 110])