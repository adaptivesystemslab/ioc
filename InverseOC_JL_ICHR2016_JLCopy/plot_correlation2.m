% plot_correlation
% h = figure('position', [   1417         169        1820         750]);
% 
% subplot(131); hold on
% colour = {'b.', 'r.', 'g.', 'm.'};
% plot(t_rmse, svd_special{length(cost_function_names)+1}(:,:), '-');
% indRange = 1:4;
% for ind = 1:length(indRange)
%     plot(t_rmse, svd_special{indRange(ind)}(:,:), colour{ind});
% end
% legend(['base' cost_function_names(indRange)]);
% % ylim([0 length(cost_function_names)]);
% 
% subplot(132); hold on
% plot(t_rmse, svd_special{length(cost_function_names)+1}(:,:), '-');
% indRange = 5:8;
% for ind = 1:length(indRange)
%     plot(t_rmse, svd_special{indRange(ind)}(:,:), colour{ind});
% end
% legend(['base' cost_function_names(indRange)]);
% % ylim([0 length(cost_function_names)]);
% 
% subplot(133); hold on
% plot(t_rmse, svd_special{length(cost_function_names)+1}(:,:), '-');
% % indRange = 9:length(cost_function_names);
% indRange = 9:length(cost_function_names);
% for ind = 1:length(indRange)
%     plot(t_rmse, svd_special{indRange(ind)}(:,:), colour{ind});
% end
% legend(['base' cost_function_names(indRange)]);
% % ylim([0 length(cost_function_names)]);
% saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankSVD_' nowStr '.fig']));
% saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankSVD_' nowStr '.png']));
% 
% close(h);

h = figure('position', [   1417         169        1820         750]);
hold on;
colour = {'b.', 'r.', 'g.', 'm.', 'k.', 'y.', 'c.', 'b.', 'r.', 'g.', 'm.', 'k.'};
plot(t_rmse, svd_special{length(cost_function_names)+1}(:,:), '-');

switch length(cost_function_names)
    case 9
        indRange = [3 8];
        
    case 12
        indRange = [3 8 11];
        
    otherwise
        indRange = 1:length(cost_function_names);
end

for ind = 1:length(indRange)
    plot(t_rmse, svd_special{indRange(ind)}(:,:) - 0.05*ind, colour{ind});
end
legend(['base' cost_function_names(indRange)]);
% ylim([0 length(cost_function_names)]);

title('SVD eigenvalues base and with items removed')
saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankSVDlarge_' nowStr '.fig']));
saveas(h, fullfile(outputPath, [currInstName '_fig22' '_rankSVDlarge_' nowStr '.png']));

close(h);

fileID = fopen(fullfile(outputPath, [currInstName '_fig22_rrefMtx.csv']), 'w+');
for ind = 1:length(rref_out)
    fmt = cell(1,length(rref_out{1}));
    fmt(:) = {',%f'};
    matrixOutputStr = horzcat(['' fmt{:} '\n']);
    fprintf(fileID, ['\n rref matrix at t= ' num2str(rref_rank(ind, 2)) ' (rank ' num2str(rref_rank(ind, 1)) ') = \n']);
    
    for ind2 = 1:length(cost_function_names)
        fprintf(fileID, '%s,', cost_function_names{ind2});
        fprintf(fileID, matrixOutputStr, rref_out{ind}(:, ind2));
    end
end

    fclose(fileID);