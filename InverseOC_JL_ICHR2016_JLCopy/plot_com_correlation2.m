% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];
ratioRMSE_plotUnique = ratioRMSE_plot(uniqueInd, :);
lengthArray = length(t_rmse);

com = feature_all.x;

upperBound = length(cost_function_names);
if length(cost_function_names) > 8
    upperBound = 8;
end

h5_2 = figure;
for ii = 5:upperBound
    for jj = 1:2
         ax(ii-4) = subplot(2, 4, (ii-5)*2 + jj); hold on
         [rmse_b, rmse_y, rmse_Rsq] = ...
             generateAndPlot2(h5_2, com(jj, end-lengthArray+1:end)', ratioRMSE_plotUnique(:, ii), [cost_function_names{ii} '_' num2str(jj)], 'b', 'x');
    end
end

% linkaxes(ax, 'xy');
