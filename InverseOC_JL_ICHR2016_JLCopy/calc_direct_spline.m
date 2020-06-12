function splineFit = calc_direct_spline(y_knot_in, param, setMode)
    % generate the knot points
%     x_traj_array = param.t([param.intermed_ind 1 1 end end]);
%     y_traj_array = [sp_ed_y zeros(size(sp_ed_y, 1), 4)];

    if nargin < 3
        setMode = param.splineType; % bformspline cubicspline polyfit ppformspline slm splinefit_lundgren piecewise_cds
    end
    
    slave_spline_mode = param.splineSlave; % slave spline for the piecewise recursive
    
    % set up the spline x, dx, ddx conditions
    [sp_ed_x, sp_ed_t, sp_ed_y, sp_ed_dx, sp_ed_dt, sp_ed_dy, sp_ed_ddx, sp_ed_ddt, sp_ed_ddy, sp_ed_dddx, sp_ed_dddt, sp_ed_dddy] ...
        = set_spline_endconditions(y_knot_in, param);
    
    switch setMode
        case {'bformspline', 'ppformspline'}
            t_traj_array = [sp_ed_t sp_ed_dt sp_ed_ddt sp_ed_dddt];
            y_traj_array = [sp_ed_y sp_ed_dy sp_ed_ddy sp_ed_dddy];

            sp_qd0 = spapi(param.spline_order, t_traj_array, y_traj_array);
            splineFit.type = 'spline';
            
            if strcmpi(setMode, 'ppformspline')
                % convert to ppform
                sp_qd0 = fn2fm(sp_qd0, 'pp');
            end

        case 'cubicspline'
            t_traj_array = sp_ed_t;
                
            if ~isempty(sp_ed_dt)
                % additional end cond are passed in
                y_traj_array = [sp_ed_dy(:, 1) sp_ed_y sp_ed_dy(:, 2)];
            else
                % no additional end conditions
                y_traj_array = sp_ed_y;
            end

            sp_qd0 = csape(t_traj_array, y_traj_array);
            splineFit.type = 'spline';
            
        case 'polyfit'
            t_traj_array = sp_ed_t;
            y_traj_array = sp_ed_y;
            
            sp_qd0 = cell(1, size(y_traj_array, 1));
            for ind_n = 1:size(y_traj_array, 1)
                sp_qd0{ind_n} = polyfit(t_traj_array, y_traj_array(ind_n, :), param.spline_order);
            end
            splineFit.type = 'polyfit';
            
        case '5th_order_poly'
            % fit for the first and last values available
            t_traj_array = sp_ed_t;
            y_traj_array = sp_ed_y;
            
            sp_qd0 = cell(1, size(y_traj_array, 1));
            for ind_n = 1:size(y_traj_array, 1)
                q0 = y_traj_array(ind_n, 1);
                qf = y_traj_array(ind_n, end);
                
                if ~isempty(sp_ed_dy)
                    dq0 = sp_ed_dy(ind_n, 1);
                    dqf = sp_ed_dy(ind_n, end);
                else
                    dq0 = 0;
                    dqf = 0;
                end
                
                if ~isempty(sp_ed_ddy)
                    ddq0 = sp_ed_ddy(ind_n, 1);
                    ddqf = sp_ed_ddy(ind_n, end);
                else
                    ddq0 = 0;
                    ddqf = 0;
                end
                
                sp_qd0{ind_n} = quintic_poly(q0, dq0, ddq0, qf, dqf, ddqf, t_traj_array);
            end
            
           splineFit.type = 'polyfit';
            
        case '7th_order_poly'

        case 'pchip'
            t_traj_array = sp_ed_t;
            y_traj_array = sp_ed_y;
            
            sp_qd0 = pchip(t_traj_array, y_traj_array);
            splineFit.type = 'spline';

        case 'piecewise_cds'
            % using this function as recursion, construct a piecewise
            % polynomial that uses each knot as the end points. enforce any
            % external dq/ddq requirements. if not explicit, then use the
            % dq/ddq from the previous stage
            
            x_traj_array = param.intermed_ind;
            y_traj_array = sp_ed_y;
            
            dy_start_prev = [];
            ddy_start_prev = [];
            dddy_start_prev = [];
            dy_end_prev = [];
            ddy_end_prev = [];
            dddy_end_prev = [];
            
            sp_qd0 = cell(1, length(x_traj_array)-1);
            for ind_knots = 1:length(x_traj_array)-1                
                % note all the constraints
                x_start_curr = x_traj_array(ind_knots);
                x_end_curr = x_traj_array(ind_knots+1);
                y_start_curr = sp_ed_y(:, ind_knots);
                y_end_curr = sp_ed_y(:, ind_knots+1);
                
                % check for dq/ddq constraints
                [dx_start_curr, dy_start_curr]     = spline_find_constraint(x_start_curr, sp_ed_dx,   sp_ed_dy,   dy_start_prev);
                [ddx_start_curr, ddy_start_curr]   = spline_find_constraint(x_start_curr, sp_ed_ddx,  sp_ed_ddy,  ddy_start_prev);
                [dddx_start_curr, dddy_start_curr] = spline_find_constraint(x_start_curr, sp_ed_dddx, sp_ed_dddy, dddy_start_prev);
                [dx_end_curr, dy_end_curr]         = spline_find_constraint(x_end_curr,   sp_ed_dx,   sp_ed_dy,   dy_end_prev);
                [ddx_end_curr, ddy_end_curr]       = spline_find_constraint(x_end_curr,   sp_ed_ddx,  sp_ed_ddy,  ddy_end_prev);
                [dddx_end_curr, dddy_end_curr]     = spline_find_constraint(x_end_curr,   sp_ed_dddx, sp_ed_dddy, dddy_end_prev);
                
                % feed into the recursive function
                param_local = param;
                param_local.intermed_ind = [x_start_curr x_end_curr];
%                 param.t_spline = (param.ti_full:(param.spline_length-1))*param.dt_spline; % just the window, what we will be calc'ing over
%                 param_local.t_spline = (x_traj_array(1):x_end_curr)*param.dt_spline;
                param_local.t_spline = param.t_spline(x_traj_array(1):x_end_curr);
                
                if param.spline_order > 2
                    param_local.sp_ed_dx = [dx_start_curr dx_end_curr];
                    param_local.sp_ed_dy = [dy_start_curr dy_end_curr];
                end
                
                if param.spline_order > 4
                    param_local.sp_ed_ddx = [ddx_start_curr ddx_end_curr];
                    param_local.sp_ed_ddy = [ddy_start_curr ddy_end_curr];
                end
                
                if param.spline_order > 6
                    param_local.sp_ed_dddx = [dddx_start_curr dddx_end_curr];
                    param_local.sp_ed_dddy = [dddy_start_curr dddy_end_curr];
                end
                
                y_knot_out = [y_start_curr y_end_curr];
                
                sp_qd0{ind_knots} = calc_direct_spline(y_knot_out, param_local, slave_spline_mode);
                [q_temp, dq_temp, ddq_temp, dddq_temp] = calc_q(sp_qd0{ind_knots}, param_local);
                
                % the end of this segment is the start of the next one
                dy_start_prev = dq_temp(:, x_end_curr);
                ddy_start_prev = ddq_temp(:, x_end_curr);
                dddy_start_prev = dddq_temp(:, x_end_curr);
            end
            splineFit.type = 'piecewise_cds';
            
            sp_qd1 = [];
            sp_qd2 = [];
            sp_qd3 = [];
    end

    switch splineFit.type
        case 'spline'
            sp_qd1 = fnder(sp_qd0,1);
            sp_qd2 = fnder(sp_qd0,2);
            sp_qd3 = fnder(sp_qd0,3);
            
        case 'polyfit'
            sp_qd1 = cell(1, size(y_traj_array, 1));
            sp_qd2 = cell(1, size(y_traj_array, 1));
            sp_qd3 = cell(1, size(y_traj_array, 1));
             for ind_n = 1:size(y_traj_array, 1)
                 sp_qd1{ind_n} = polyder_num(sp_qd0{ind_n});
                 sp_qd2{ind_n} = polyder_num(sp_qd1{ind_n});
                 sp_qd3{ind_n} = polyder_num(sp_qd2{ind_n});
             end
             
        case 'piecewise_cds'
            % no need, since the differentiation for each segment is
            % already completed in the recursion
                
    end
 
    splineFit.x_knot = param.intermed_ind;
    splineFit.t_spline = param.t_spline;
    splineFit.sp_ed_y = sp_ed_y;
    splineFit.sp_qd0 = sp_qd0;
    splineFit.sp_qd1 = sp_qd1;
    splineFit.sp_qd2 = sp_qd2;
    splineFit.sp_qd3 = sp_qd3;
    
    if 0
        q =     fnval(splineFit.sp_qd0, param.t_spline);     figure; plot(q');
        dq_sp =     fnval(splineFit.sp_qd1, param.t_spline);
        ddq_sp =     fnval(splineFit.sp_qd2, param.t_spline);
        
        dq_diff = calcDeriv(q, param.dt_spline);
        ddq_diff = calcDeriv(dq, param.dt_spline);
        
        q =     fnval(sp_qd0{1}, param.t_spline);     figure; plot(q');
        
         [q, dq, ddq, dddq] = calc_q(splineFit, param);
    end
end

function [dx_curr, dy_curr]= spline_find_constraint(x_curr, sp_ed_dx, sp_ed_dy, dy_prev)
    ind_find = find(sp_ed_dx == x_curr);
    
    if ~isempty(ind_find)
        % found a constraint passed in from outside
        dx_curr = x_curr;
        dy_curr = sp_ed_dy(:, ind_find);
    elseif ~isempty(dy_prev)
        % use the previous constraint
        dx_curr = x_curr; 
        dy_curr = dy_prev;
    else
        dx_curr = [];
        dy_curr = [];
    end
end

function [sp_ed_x, sp_ed_t, sp_ed_y, sp_ed_dx, sp_ed_dt, sp_ed_dy, sp_ed_ddx, sp_ed_ddt, sp_ed_ddy, sp_ed_dddx, sp_ed_dddt, sp_ed_dddy] ...
    = set_spline_endconditions(y_knot_in, param)
    
    if iscell(y_knot_in)
        sp_ed_x = param.intermed_ind;
        sp_ed_t = param.t_spline(param.intermed_ind);
        sp_ed_y = y_knot_in{1};
    else
        sp_ed_x = param.intermed_ind;
        sp_ed_t = param.t_spline(param.intermed_ind);
        sp_ed_y = y_knot_in;
    end

    % additional end cond are passed in
    if iscell(y_knot_in)
        sp_ed_dx = param.intermed_ind;
        sp_ed_dt = param.t_spline(param.intermed_ind);
        sp_ed_dy = y_knot_in{2};
    elseif isfield(param, 'sp_ed_dx') && ~isempty(param.sp_ed_dx) && param.spline_order > 2
        sp_ed_dx = param.sp_ed_dx;
        sp_ed_dt = param.t_spline(param.sp_ed_dx);
        sp_ed_dy = param.sp_ed_dy;
    else
        sp_ed_dx = [];
        sp_ed_dt  = [];
        sp_ed_dy = [];
    end
    
    if iscell(y_knot_in)
        sp_ed_ddx = param.intermed_ind;
        sp_ed_ddt = param.t_spline(param.intermed_ind);
        sp_ed_ddy = y_knot_in{3};
    elseif isfield(param, 'sp_ed_ddx') && ~isempty(param.sp_ed_ddx) && param.spline_order > 4
        sp_ed_ddx = param.sp_ed_ddx;
        sp_ed_ddt = param.t_spline(param.sp_ed_ddx);
        sp_ed_ddy = param.sp_ed_ddy;
    else
        sp_ed_ddx = [];
        sp_ed_ddt = [];
        sp_ed_ddy = [];
    end
    
    if isfield(param, 'sp_ed_dddx') && ~isempty(param.sp_ed_dddx) && param.spline_order > 6
        sp_ed_dddx = param.sp_ed_dddx;
        sp_ed_dddt = param.t_spline(param.sp_ed_dddx);
        sp_ed_dddy = param.sp_ed_dddy;
    else
        sp_ed_dddx = [];
        sp_ed_dddt = [];
        sp_ed_dddy = [];
    end
end