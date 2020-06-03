function h = plot_direct(param, y_opt_sept, dy_opt_sept, ddy_opt_sept, q, dq, ddq, tau, dtau, grf)
    h = figure('position', [   1415         230        1711         840]);
    subplot(231);
    plot(param.t, y_opt_sept', 'x'); hold on
    plot(param.t, q', 'o'); 
    plot(param.t, q', '-');
%     line(param.t([ind ind]), ylim);
%     line(param.t([1 1]), ylim);
%     line(param.t([length(y_opt_sept) length(y_opt_sept)]), ylim);
    for i = 1:length(param.x_knot)
        line([param.x_knot(i) param.x_knot(i)], ylim, 'Color', 'b');
    end
    for i = 1:length(param.const_x)
        line(param.t([param.const_x(i) param.const_x(i)]), ylim, 'Color', 'r');
    end
    title('q');

    subplot(232);
    plot(param.t, dy_opt_sept', 'x'); hold on
    plot(param.t,dq', 'o'); 
    plot(param.t,dq', '-');
    title('dq');    
    
    subplot(233);
    plot(param.t, ddy_opt_sept', 'x'); hold on
    plot(param.t,ddq', 'o'); 
    plot(param.t,ddq', '-');
    title('ddq');
    
    subplot(234);
    plot(param.t, tau');
    title('torque');
    
    subplot(235);
    plot(param.t, dtau');
    title('diff torque');
    
    subplot(236);
    plot(param.t, grf');
    title('grf');
end