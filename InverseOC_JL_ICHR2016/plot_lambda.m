% % % % plot data
% % % % [~, uniqueInd] = unique(t_rmse);
% % % ax = [];
% % % 
% % % if exist('h14l', 'var') && h14l > 0
% % %     figure(h14l);
% % % else
% % %     h14l = figure;
% % % end
% % % 
% % % ax(1) = subplot(2, 1, 1); hold on
% % % plot(feature_full.t, q_opt_plot_merge);
% % % plot(t_recon_plot_merge, q_recon_plot_merge, '.', 'LineWidth', 1);
% % % title(['Trajectory and Reconstruction (RMSE (blended, windowed) : ' ...
% % %     num2str(rmse_report.blended_rmse_belowThreshold) ', ' num2str(rmse_report.windowed_rmse_belowThreshold_mean) ' \pm ' num2str(rmse_report.windowed_rmse_belowThreshold_std) ')']);
% % % ylabel('Joint angle [rad]','Interpreter','Latex');
% % % plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0)), ...
% % %     zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 0))) - 2, 'rx');
% % % plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1)), ...
% % %     zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass == 1))) - 2, 'yx');
% % % plot(feature_full.t(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1)), ...
% % %     zeros(size(find(resnromAll_lsqlin_const_minRMSE_array_belowThresPass > 1))) - 2, 'gx');
% % % plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 1.05, -2.05);
% % % 
% % % % switch run_mode
% % % %     case 'win'
% % % %         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 2, -2);
% % % % end
% % % 
% % % hold on
% % % headers = {'Ankle $$q_{obs}$$', 'Knee $$q_{obs}$$', 'Hip $$q_{obs}$$', 'Ankle $$\hat{q}_{obs}$$', 'Knee $$\hat{q}_{obs}$$', 'Hip $$\hat{q}_{obs}$$', 'Pass = 0', 'Pass = 1', 'Pass \> 1'};
% % % % legend(headers,'Interpreter','Latex');
% % % 
% % % 
% % % for ind_lambdaCounter = 1:length(lambda_name)
% % %     ax(2) = subplot(2, 1, 2); 
% % %     % title('Recovered Contribution Weights');
% % %     % indToUse = [1:3 10:12 19:21];
% % %     indToUse = ind_lambdaCounter;
% % %     % lambda_array_to_use = lambda_array(:, indToUse);
% % %     lambda_array_to_use = avgLambdaArray_belowThres(:, indToUse);
% % % %     norm_lambda = normVector(lambda_array_to_use);
% % %     area(feature_full.t, lambda_array_to_use); 
% % % %     plot(feature_full.t, norm_lambda, 'k.'); 
% % %     shading flat
% % %     % ylim([0 1])
% % %     plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, 0.45, -0.45);
% % %     % ylabel('Recovered $$\hat{c}$$','Interpreter','Latex');
% % %     xlabel('Time [s]');
% % %     title(lambda_name{ind_lambdaCounter});
% % %     % legend(cost_function_names_sorted);
% % %     ylim([-0.5 0.5]);
% % %     
% % %     linkaxes(ax, 'x');
% % %     
% % %     saveas(h14l, fullfile(outputPath, [currInstName '_fig15_lambda' num2str(ind_lambdaCounter) '_' nowStr '.fig']));
% % %     saveas(h14l, fullfile(outputPath, [currInstName '_fig15_lambda' num2str(ind_lambdaCounter) '_' nowStr '.png']));
% % % end
% % % 
% % % 
% % % % xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);

% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h14l', 'var') && h14l > 0
    figure(h14l);
else
    h14l = figure;
end

for ind_lambdaCounter = 1:3:length(lambda_name)
    ax(1) = subplot(3, 1, 1); 
    indToUse = ind_lambdaCounter;
    lambda_array_to_use = avgLambdaArray_belowThres(:, indToUse);
    area(feature_full.t, lambda_array_to_use); 
    shading flat
    plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
    xlabel('Time [s]');
    title(lambda_name{indToUse});
%     ylim([-0.5 0.5]);
    
    ax(2) = subplot(3, 1, 2); 
    indToUse = ind_lambdaCounter+1;
    lambda_array_to_use = avgLambdaArray_belowThres(:, indToUse);
    area(feature_full.t, lambda_array_to_use); 
    shading flat
    plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
    xlabel('Time [s]');
    title(lambda_name{indToUse});
%     ylim([-0.5 0.5]);  
    
    ax(3) = subplot(3, 1, 3); 
    indToUse = ind_lambdaCounter+2;
    lambda_array_to_use = avgLambdaArray_belowThres(:, indToUse);
    area(feature_full.t, lambda_array_to_use); 
    shading flat
    plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
    xlabel('Time [s]');
    title(lambda_name{indToUse});
%     ylim([-0.5 0.5]);      
    
    linkaxes(ax, 'xy');
    
    saveas(h14l, fullfile(outputPath, [currInstName '_fig15_lambdaSet' num2str(ind_lambdaCounter) '_' nowStr '.fig']));
    saveas(h14l, fullfile(outputPath, [currInstName '_fig15_lambdaSet' num2str(ind_lambdaCounter) '_' nowStr '.png']));
end


% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);
