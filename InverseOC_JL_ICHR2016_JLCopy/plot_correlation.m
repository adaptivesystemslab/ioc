% plot_correlation
for ind = 1:length(cost_function_names)
    h1 = figure('position', [   1417         169        1820         750]);
    
    subplot(131);
    plot(t_rmse, correlationArray_special{ind}(:, 1:4));
    title([num2str(ind) ': ' cost_function_names{ind}]);
    legend(cost_function_names(1:4));
    
    subplot(132);
    plot(t_rmse, correlationArray_special{ind}(:, 5:8));
    legend(cost_function_names(5:8));
    
    subplot(133);
    plot(t_rmse, correlationArray_special{ind}(:, 9:end));
    legend(cost_function_names(9:end));
    
    saveas(h1, fullfile(outputPath, [currInstName '_fig21_' num2str(ind) '_corr_' nowStr '.fig']));
    saveas(h1, fullfile(outputPath, [currInstName '_fig21_' num2str(ind) '_corr_' nowStr '.png']));
    
    close(h1);
end



% plot_correlation
h = figure('position', [   1417         169        1820         750]);

subplot(131); hold on
colour = {'b', 'r', 'g', 'm'};
plot(t_rmse, rankGrad_special{length(cost_function_names)+1}(:,:), '.');
indRange = 1:4;
for ind = 1:length(indRange)
    plot(t_rmse, rankGrad_special{indRange(ind)}(:,:), colour{ind});
end
legend(cost_function_names(indRange));
% ylim([0 length(cost_function_names)]);

subplot(132); hold on
plot(t_rmse, rankGrad_special{length(cost_function_names)+1}(:,:), '.');
indRange = 5:8;
for ind = 1:length(indRange)
    plot(t_rmse, rankGrad_special{indRange(ind)}(:,:), colour{ind});
end
legend(cost_function_names(indRange));
% ylim([0 length(cost_function_names)]);

subplot(133); hold on
plot(t_rmse, rankGrad_special{length(cost_function_names)+1}(:,:), '.');
indRange = 9:length(cost_function_names);
for ind = 1:length(indRange)
    plot(t_rmse, rankGrad_special{indRange(ind)}(:,:), colour{ind});
end

legend(cost_function_names(indRange));
% ylim([0 length(cost_function_names)]);
saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankG_' nowStr '.fig']));
saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankG_' nowStr '.png']));

close(h);



% plot_correlation
h = figure('position', [   1417         169        1820         750]);

subplot(131); hold on
colour = {'b', 'r', 'g', 'm'};
indRange = 1:4;
for ind = 1:length(indRange)
    plot(t_rmse, rankA_special{indRange(ind)}(:,:), colour{ind});
end
legend(cost_function_names(indRange));
% ylim([0 length(cost_function_names)]);

subplot(132); hold on
indRange = 5:8;
for ind = 1:length(indRange)
    plot(t_rmse, rankA_special{indRange(ind)}(:,:), colour{ind});
end
legend(cost_function_names(indRange));
% ylim([0 length(cost_function_names)]);

subplot(133); hold on
indRange = 9:length(cost_function_names);
for ind = 1:length(indRange)
    plot(t_rmse, rankA_special{indRange(ind)}(:,:), colour{ind});
end
legend(cost_function_names(indRange));
% ylim([0 length(cost_function_names)]);

saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankA_' nowStr '.fig']));
saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankA_' nowStr '.png']));

close(h);
