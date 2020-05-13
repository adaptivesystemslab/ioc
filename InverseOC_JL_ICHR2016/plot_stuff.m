% plot more stuff



for ind_plotCount = 1:ccostCount
    h8(ind_plotCount) = figure;
    ax(1) = subplot(3, 1, 1); hold on
    plot(feature_full.t, q_opt_plot_merge);
    plot(t_recon_plot_merge, q_recon_plot_merge, '.', 'LineWidth', 1);
    title(['RMSE Mean: ' num2str(rmse_report.meanRMSE) ', STD: ' num2str(rmse_report.stdRMSE), ', Range: ' num2str(rmse_report.rangeRMSE)]);
    
    switch run_mode
        case 'sim'
            plotBoxes(gcf, feature_full.t(indToUse_window(:, 1)), feature_full.t(indToUse_window(:, 2)));
            
        case 'win'
            plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
    end
    
    ax(2) = subplot(3, 1, 2); hold on
    title(['Contribution weight for ' cost_function_names{ind_plotCount}]);
    plot(t_rmse(pts_underConstThres), zeros(size(find(pts_underConstThres))) , 'gx');
    area(t_rmse(uniqueInd), cAll_plot(uniqueInd, ind_plotCount));
   
%     plot(t_rmse(pts_aboveConstThres), zeros(size(find(pts_aboveConstThres)))  , 'rx');

    
%     plot(t_rmse(pts_aboveDiffThres), zeros(size(find(pts_aboveDiffThres)))  -0.0003, 'rx');
%     plot(t_rmse(pts_underDiffThres), zeros(size(find(pts_underDiffThres)))  -0.0003, 'gx');
    
    ax(3) = subplot(3, 1, 3); hold on
    title(['Mean window weight for ' cost_function_names{ind_plotCount}]);
    plot(feature_full.t(pts_aboveConstThres_averageWindow), zeros(size(find(pts_aboveConstThres_averageWindow))) , 'gx');
    area(feature_full.t(:), avgWeightArray(:, ind_plotCount));
    
%     plot(t_rmse(pts_aboveConstThres), zeros(size(find(pts_aboveConstThres)))  , 'rx');

    
%     plot(t_rmse(pts_aboveDiffThres), zeros(size(find(pts_aboveDiffThres)))  -0.0003, 'rx');
%     plot(t_rmse(pts_underDiffThres), zeros(size(find(pts_underDiffThres)))  -0.0003, 'gx');
    
    linkaxes(ax, 'x');
    % xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);
end

