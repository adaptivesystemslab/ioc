% plot data
% [~, uniqueInd] = unique(t_rmse);
ax = [];

if exist('h11_2', 'var') && h11_2 > 0
    figure(h11_2);
else
    h11_2 = figure;
end

ax(1) = subplot(2, 1, 1); hold on
% title('Averaged Normalized Cost Weights');
bar(feature_full.t, avgWeightArray_belowThres, 'stacked'); 
shading flat
legend(cost_function_names);
% xlabel('Time [s]');
ylabel('Normalized Cost Weights');

ylim([-0.05 1.05])

ax(2) = subplot(2, 1, 2); hold on
area(feature_full.t, resnromAll_lsqlin_const_minRMSE_array_belowThres(:, :), 'FaceColor', [0 0 1]); 
la = ylim;

title(['Averaged Residual Norm:' num2str(mean(resnromAll_lsqlin_const_minRMSE_array_belowThres))]);
plot(xlim, [constThres constThres], 'k');
xlabel('Time [s]');
ylabel('Averaged Residual Norm');
ylim(la);
% upperbound = min([constThres*5 max(resnromAll_lsqlin_const_minRMSE_array_belowThres)]);

% ylim([0 constThres*5]);

linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);