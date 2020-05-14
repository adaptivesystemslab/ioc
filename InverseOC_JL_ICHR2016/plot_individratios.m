% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

h3_1 = figure;
upperBound = 4;

if size(ratioRMSE_plot, 2) < upperBound
    upperBound = size(ratioRMSE_plot, 2);
end

for ii = 1:upperBound
    ax(ii) = subplot(2, 2, ii); hold on
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

linkaxes(ax, 'xy');
xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);
ylim([-10 110])