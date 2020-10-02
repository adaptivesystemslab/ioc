% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];
ratioRMSE_plotUnique = ratioRMSE_plot(uniqueInd, :);
lengthArray = length(t_rmse);

% com = sqrt(sum(feature_all.com.^2));
com = feature_all.com;

upperBound = 4;
if size(ratioRMSE_plot, 2) < upperBound
    upperBound = size(ratioRMSE_plot, 2);
end

h5_1 = figure;
for ii = 1:upperBound
    for jj = 1:2
         ax(ii) = subplot(2, 4, (ii-1)*2 + jj); hold on
         [rmse_b, rmse_y, rmse_Rsq] = ...
             generateAndPlot2(h5_1, com(jj, end-lengthArray+1:end)', ratioRMSE_plotUnique(:, ii), [cost_function_names{ii} '_' num2str(jj)], 'b', 'x');
    end
end

% linkaxes(ax, 'xy');
