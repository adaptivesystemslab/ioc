function [rmse_b, rmse_y, rmse_Rsq] = generateAndPlot(h, rmse_ioc, c_ioc_ratio_error, minRmseInd, minResnormInd, name)
    [rmse_b, rmse_y, rmse_Rsq] = calcRegression(rmse_ioc, c_ioc_ratio_error);

    hold on
    plot(rmse_ioc, c_ioc_ratio_error, 'bx', 'LineWidth', 5, 'DisplayName', 'All');
    plot(rmse_ioc(minRmseInd), c_ioc_ratio_error(minRmseInd), 'ro', 'LineWidth', 5, 'DisplayName', 'rmse');
    plot(rmse_ioc(minResnormInd), c_ioc_ratio_error(minResnormInd), 'go', 'LineWidth', 5, 'DisplayName', 'resnorm');
    plot(rmse_ioc, rmse_y);
    title([name ', R=' num2str(rmse_Rsq)]);
    ylim([-1 201])
end