function [feature_use, J_array, J_contrib, output_direct] = main_direct(c_cost, c_const, feature_win, param, param_opt)
    % direct problem - given a cost function, generate trajectory. the
    % states is the knot points
    
    displaySetting = 'iter-detailed'; % 'final' 'notify' 'off' 'iter-detailed'
    
    param = update_intermed_ind(param, param_opt.knot_count);
    
    % using a polynomial least squares to roughly generate an initial cond.
    % the function is a 7th order fit, but assumes dq=ddq=dddq=0
    [y_init, dy_init, ddy_init] = generateTrajectory(param.const_y, param.const_x, param.t_spline, zeros(size(param.const_y)), zeros(size(param.const_y)), param);
    
    const_y_delta = (180/pi) * (param.const_y(:, 2) - param.const_y(:, 1));
    
    %    y_opt_sept = calc_trajectory_septic(param.t, [param.const_y{1} param.const_y{2}]);
    y_opt_base{1} = y_init(:, param.intermed_ind); % pull out the knot points
    y_opt_base{2} = dy_init(:, param.intermed_ind); % dy
    y_opt_base{3} = ddy_init(:, param.intermed_ind); % ddy

    [y_opt_vec, param.matVec_struct] = mat2vec_multi(y_opt_base); % converting the state matrix (easier to read) to a vector (fmincon req)
 
    splineFit.q = calc_direct_spline(y_opt_base, param);
    feature_use = calc_features(splineFit, feature_win, param);
    [J, J_contrib, J_array] = calc_direct_cost(c_cost, feature_use, param);
    [ceq, ceq_x, ceq_dx, ceq_ddx] = calc_direct_const(c_const, feature_use, param);
     
    % fmincon setup
    cost_func_handle =  @(y_opt) direct_cost_func(y_opt, c_cost, c_const, feature_win, param);
    const_func_handle = @(y_opt) direct_const_func(y_opt, c_cost, c_const, feature_win, param);
    
    cost_funcZero_handle = @(y_opt) 0;
    
    Mesoptions = optimset('Jacobian','on','Algorithm','interior-point','FunValCheck','on',...
        'MaxFunEvals',1e8,'MaxIter',param_opt.maxiter, 'Display',displaySetting, ... % iter-detailed
        'TolFun',param_opt.TolFun,'TolX',param_opt.TolX,'TolCon',param_opt.TolCon);
    
%     % may need to bootstrap a solution
%     [y_opt_vec, cost_J, exitflag] = fmincon(cost_funcZero_handle, y_opt_vec,[],[],[],[],[],[],const_func_handle,Mesoptions);
    
    % actual solver
    [y_out_vec, cost_J, exitflag, output, lambda, grad, hessian] = fmincon(cost_func_handle, y_opt_vec,[],[],[],[],[],[],const_func_handle,Mesoptions);
    exitflag
    
    % calculate the results to pass out of the function
    y_opt_mat = vec2mat_multi(y_out_vec, param.matVec_struct);   
    
    splineFit.q = calc_direct_spline(y_opt_mat, param);
    feature_use = calc_features(splineFit, feature_win, param);
    [J, J_contrib, J_array, J_debug] = calc_direct_cost(c_cost, feature_use, param);
    [ceq, ceq_x, ceq_dx, ceq_ddx] = calc_direct_const(c_const, feature_use, param);
    
    feature_use.xknot_ind = param.intermed_ind;
    feature_use.xknot_t = param.x_knot;
    
    if 0
        dq = calcDeriv(feature_use.q, param.dt_spline);
        ddq = calcDeriv(dq, param.dt_spline);
%         dddq = calcDeriv(ddq, param.dt_spline);
        
        h = figure;
        subplot(211);
        plot(feature_use.q');
        hold on
        plot(param.const_x, param.const_y', 'x',  'MarkerSize', 10);
        title(['norm(violation of ceq-x) = ' num2str(norm(ceq_x))]);
        ylabel('q');
        xlabel(['exitflag = ' num2str(exitflag)]);
        
        subplot(212);
        plot(ceq_ddx');
        ylabel('ceq-ddx');
        title(['norm(violation of ceq-ddx) = ' num2str(norm(ceq_ddx))]);
        saveas(h, ['saved\' datestr(now, 'yyyymmddHHMMSS')]);
        
%         xlabel(['exitflag = ' num2str(exitflag)]);
        
%         figure;
%         subplot(221);
%         plot(feature_use.q'); title('q');
%         subplot(222);
%         plot(feature_use.dq'); title('dq');
%         hold on; plot(dq', 'r');
%         subplot(223);
%         plot(feature_use.ddq'); title('ddq');
%         hold on; plot(ddq', 'r');
%         subplot(224);
%         plot(feature_use.dddq'); title('dddq');
%         hold on; plot(dddq', 'r');
    end
   

    output_direct.grad = grad;
    output_direct.exitflag = exitflag;
    output_direct.J_array = J_array;
    output_direct.J_contrib = J_contrib;
    output_direct.c = 0;
    output_direct.J = J;
    output_direct.J_breakdown = J_debug;
    output_direct.maxceq = max(abs(ceq));
    output_direct.ceq_x = norm(ceq_x);
    output_direct.ceq_dx = norm(ceq_dx);
    output_direct.ceq_ddx = norm(ceq_ddx);
    output_direct.ceq_x_array = (ceq_x);
    output_direct.ceq_dx_array = (ceq_dx);
    output_direct.ceq_ddx_array = (ceq_ddx);    
end

% cost function
function J = direct_cost_func(y_opt_vec, c_cost, c_const, feature_win, param)
    y_opt_mat = vec2mat_multi(y_opt_vec, param.matVec_struct); % convert the vector back into matrix form for reading
    splineFit.q = calc_direct_spline(y_opt_mat, param); % calc the spline from knot points (ie the variables to opt for)
    
    feature_use = calc_features(splineFit, feature_win, param);

    % cost functions
    J = calc_direct_cost(c_cost, feature_use, param);
end

% constrant function
function [c, ceq] = direct_const_func(y_opt_vec, c_cost, c_const, feature_win, param)
    y_opt_mat = vec2mat_multi(y_opt_vec, param.matVec_struct);
    splineFit.q = calc_direct_spline(y_opt_mat, param);
    
    feature_use = calc_features(splineFit, feature_win, param);
    
    % no inequality constraints
    c = []; 

    % equality constrains are postural
    ceq = calc_direct_const(c_const, feature_use, param);
end

