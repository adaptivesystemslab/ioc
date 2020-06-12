function [rmse_b, rmse_y, rmse_Rsq] = generateAndPlot2(h, x, y, name, color, sym)
    [rmse_b, rmse_y, rmse_Rsq] = calcRegression(x, y);

    hold on
    plot(x, y, [color sym], 'LineWidth', 5, 'DisplayName', 'All');
%     plot(rmse_ioc(minRmseInd), c_ioc_ratio_error(minRmseInd), 'ro', 'LineWidth', 5, 'DisplayName', 'rmse');
%     plot(rmse_ioc(minResnormInd), c_ioc_ratio_error(minResnormInd), 'go', 'LineWidth', 5, 'DisplayName', 'resnorm');
    plot(x, rmse_y, color);
    title([name ', R=' num2str(rmse_Rsq)]);
%     ylim([-1 201])
end