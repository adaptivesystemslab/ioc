% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h11_3', 'var') && h11_3 > 0
    figure(h11_3);
else
    h11_3 = figure;
end

ax(1) = subplot(3, 1, 1); hold on
% title('Averaged Normalized Cost Weights');
bar(feature_full.t, avgWeightArray_belowThres, 'stacked'); 
shading flat
legend(cost_function_names);
ylim([-0.05 1.05])
ylabel('Basis Weights');

ax(2) = subplot(3, 1, 2); hold on
area(feature_full.t, resnromAll_lsqlin_const_minRMSE_array(:, :), 'FaceColor', [0 0 1]); 

% title('Averaged Residual Norm');
plot(xlim, [constThres constThres], 'k');
title(['mu = ' num2str(mean(resnromAll_lsqlin_const_minRMSE_array)) ', std = ' num2str(std(resnromAll_lsqlin_const_minRMSE_array))]);
% xlabel('Time [s]');
ylabel('Pre-thres Residual');

ax(3) = subplot(3, 1, 3); hold on
area(feature_full.t, resnromAll_lsqlin_const_minRMSE_array_belowThres(:, :), 'FaceColor', [0 0 1]); 
title(['mu = ' num2str(mean(resnromAll_lsqlin_const_minRMSE_array_belowThres)) ', std = ' num2str(std(resnromAll_lsqlin_const_minRMSE_array_belowThres))]);

% title('Averaged Residual Norm');
plot(xlim, [constThres constThres], 'k');
xlabel('Time [s]');
ylabel('Post-thres Residual');



linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);