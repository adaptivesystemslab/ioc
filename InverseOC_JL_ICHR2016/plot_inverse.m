figure; plot(param.t_spline, feature_opt.q', 'o', 'DisplayName', 'Original');
    hold on
    line([param.x_knot' param.x_knot'], ylim, 'Color', 'k');
    % plot(param.x_knot', param.y_knot', 'o', 'DisplayName', 'Knot points');
    for ii = 1:size(param.spline_fit, 1)
        for jj = 1:size(param.spline_fit, 2)
            toPlot = fnval(param.spline_fit(ii, jj).sp_qd0,param.t);
            plot(param.t', toPlot, '-', 'DisplayName', 'SplineRecon');
        end
    end
    title('Spline - q d0');
%     legend show