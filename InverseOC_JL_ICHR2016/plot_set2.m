% plot data
h = figure;
ax(1) = subplot(4, 4, [1 2 3 5 6 7]); hold on
plot(t_full, q_opt_plot_merge);
plot(t_recon_plot, q_recon_plot_merge, '.', 'LineWidth', 1);
title(['RMSE Mean: ' num2str(meanRMSE) ', STD: ' num2str(stdRMSE), ', Range: ' num2str(rangeRMSE)]);
        
switch mode
    case 'sim'        
        plotBoxes(gcf, t_full(indToUse_window(:, 1)), t_full(indToUse_window(:, 2)));
        
        ax(5) = subplot(4, 4, [9 13]); hold on
        plot(t_full, J_full(:, 1), 'Color', [0.6 0.6 0.1], 'LineWidth', 5, 'DisplayName', 'J1 GT');
        plot(t_rmse, ratioRMSE_plot(:, 1), 'Color', [0.2 0.2 0.1], 'LineWidth', 2, 'DisplayName', 'J1 Reg');
        ylim([-10 110]);
        title('J ratio');
        
        ax(6) = subplot(4, 4, [10 14]); hold on
        plot(t_full, J_full(:, 2), 'Color', [0.6 0.1 0.6], 'LineWidth', 5, 'DisplayName', 'J2 GT');
        plot(t_rmse, ratioRMSE_plot(:, 2), 'Color', [0.2 0.1 0.2], 'LineWidth', 2, 'DisplayName', 'J2 Reg');
        ylim([-10 110]);
        title('J2 ratio');
        
        ax(7) = subplot(4, 4, [11 15]); hold on
        plot(t_full, J_full(:, 3), 'Color', [0.1 0.6 0.6], 'LineWidth', 5, 'DisplayName', 'J3 GT');
        plot(t_rmse, ratioRMSE_plot(:, 3), 'Color', [0.1 0.2 0.2], 'LineWidth', 2, 'DisplayName', 'J3 Reg');
        ylim([-10 110]);
        title('J3 ratio');
        
    case 'win'
        plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
        
        ax(5) = subplot(4, 4, [9 10 11 13 14 15]); hold on

        if 0
            plot(t_rmse, ratioRMSE_plot(:, 1), 'rx', 'LineWidth', 2, 'DisplayName', '1) ddq');
            plot(t_rmse, ratioRMSE_plot(:, 2), 'bo', 'LineWidth', 2, 'DisplayName', '2) ddx');
            plot(t_rmse, ratioRMSE_plot(:, 3), 'gx', 'LineWidth', 2, 'DisplayName', '3) tau');
            
            plot(t_rmse, ratioRMSE_plot(:, 1), 'r');
            plot(t_rmse, ratioRMSE_plot(:, 2), 'b');
            plot(t_rmse, ratioRMSE_plot(:, 3), 'g');
            legend show
        else
%             ratioRMSE_set(:, 1) = ratioRMSE_plot(:, 1);
%             ratioRMSE_set(:, 2) = ratioRMSE_set(:, 1) + ratioRMSE_plot(:, 2);
%             ratioRMSE_set(:, 3) = ratioRMSE_set(:, 2) + ratioRMSE_plot(:, 3);
%             
%             plot(t_rmse, ratioRMSE_set(:, 1), 'r');
%             plot(t_rmse, ratioRMSE_set(:, 2), 'b');
%             plot(t_rmse, ratioRMSE_set(:, 3), 'g');

            title('Contribution ratio');
            [~, uniqueInd] = unique(t_rmse);
            bar(t_rmse(uniqueInd), ratioRMSE_plot(uniqueInd, :), 'stacked')
            shading flat
            legend('ddq', 'ddx', 'tau');
        end
%          = plotBoxes(h, currTimeStart, currTimeEnd, colorToUse, offset, maxY, minY)
        plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', [0 0], 105, -5);
        
        ylim([-10 110])
end

ax(2) = subplot(4, 4, [4 8]); hold on
plot(t_rmse, rmseMin_plot, 'k', 'LineWidth', 2); ylimGet = ylim; % plot first time to generate the ylim data
bar(t_rmse, rmseAll_plot, 'BarWidth', 1, 'BarLayout', 'grouped'); ylim(ylimGet); shading flat
plot(t_rmse, rmseMin_plot, 'k', 'LineWidth', 2);
title('rmse');

ax(3) = subplot(4, 4, 12); hold on
plot(t_rmse, resnormMin_plot, 'k', 'LineWidth', 2); ylimGet = ylim;
bar(t_rmse, resnormAll_plot, 'BarWidth', 1, 'BarLayout', 'grouped'); ylim(ylimGet); shading flat
plot(t_rmse, resnormMin_plot, 'k', 'LineWidth', 2);
title('resnorm');

ax(4) = subplot(4, 4, 16); hold on
plot(t_rmse, condMin_plot, 'k', 'LineWidth', 2); ylimGet = ylim;
bar(t_rmse, condAll_plot, 'BarWidth', 1, 'BarLayout', 'grouped'); ylim(ylimGet); shading flat
plot(t_rmse, condMin_plot, 'k', 'LineWidth', 2);
title('cond');

linkaxes(ax, 'x');