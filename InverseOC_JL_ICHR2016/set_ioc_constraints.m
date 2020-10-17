function [param, set_ioc_check_struct] = set_ioc_constraints(feature_win, feature_full, startTime, param, constMode)

if nargin == 3
    constMode = 'ioc';
end

switch constMode
    case 'ioc'
        cf_q_constraint =    param.ioc_cf_q_constraint;
        cf_q_constraint_y =  param.ioc_cf_q_constraint_y;
        cf_dq_constraint =   param.ioc_cf_dq_constraint;
        cf_dq_constraint_y = param.ioc_cf_dq_constraint_y;

    case 'doc_rmse'
        cf_q_constraint =    param.doc_rmse_cf_q_constraint;
        cf_q_constraint_y =  param.doc_rmse_cf_q_constraint_y;
        cf_dq_constraint =   param.doc_rmse_cf_dq_constraint;
        cf_dq_constraint_y = param.doc_rmse_cf_dq_constraint_y;
end

knot_locations = param.ioc_knot_locations;
add_knot_as_cf_constraints = param.ioc_add_knot_as_cf_constraints;
sp_dq_endcond = param.ioc_sp_dq_endcond;

t_win_findMatchArray = feature_win.t;
q_win_findMatchArray = 1:length(feature_win.t);

% determine the knot points
switch knot_locations
    case 'balanced'
        length_traj = size(param.t_spline, 2);
        skip_count = ceil(length_traj/param.ioc.knot_count);
        intermed_ind = unique([1:skip_count:length_traj-1 length_traj]); % this ensures there is a knot point at the end, and on constraints

    case 'variable'  
        length_traj = size(param.t_spline, 2);
        t_full = feature_full.t_knot;
        q_full = feature_full.t_knot;
        intermed_ind = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        intermed_ind = unique([1 intermed_ind length_traj]);
        
        % if knot points are too close together, remove them
        if intermed_ind(2) - intermed_ind(1) < 5
            intermed_ind(2) = 0;
        end
        
        if intermed_ind(end) - intermed_ind(end-1) < 5
            intermed_ind(end-1) = 0;
        end
        
        intermed_ind = intermed_ind(intermed_ind > 0);
        
    case 'double'
        length_traj = size(param.t_spline, 2);
        skip_count = ceil(length_traj/(param.ioc.knot_count*2));
        intermed_ind = unique([1:skip_count:length_traj-1 length_traj]); % this ensures there is a knot point at the end, and on constraints
        
    case 'matchRand'
        t_full = feature_full.t_knot;
        q_full = feature_full.t_knot;
        intermed_ind = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        
        % randomize the inner points a bit
        intermed_ind(2:end) = intermed_ind(2:end) + 1;
end

switch cf_q_constraint
    case 'fixed'
        %             const_y{1} = feature_win.q(:, param.const_x(1)).*repmat(param.coeff_q,     1, 1); % the constraints are passed to the IOC and is compared against the normalized q, so should be normalized first
        %             const_y{2} = feature_win.q(:, param.const_x(2)).*repmat(param.coeff_q,     1, 1);
        %             const_y{3} = feature_win.q(:, param.const_x(3)).*repmat(param.coeff_q,     1, 1);
        const_x = param.const_x;
        
    case 'variable'
        % find all the const_x in the window
        t_target = feature_full.t_cf_q_const;
        q_target = 1:length(t_target);
        const_x = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_target, q_target, param);
        
    case {'startend', 'edge'}
        const_x = [1 param.spline_length];
        
    case {'start'}
        const_x = [1];
        
    case {'end'}
        const_x = [param.spline_length];
        
    case {'startmid1end', 'startmidend'}
        const_x = [1 ceil(param.spline_length/2) param.spline_length];
        
    case {'startmid2end', 'startmidmidend'}
        const_x = [1 ceil(param.spline_length/3) ceil(2*param.spline_length/3) param.spline_length];
        
    case 'startmid1endonknot'
        mid1 = ceil(param.spline_length/2)+1;
        [findVal1, findInd1] = findClosestValue(mid1, intermed_ind, 'above'); % find the closest knot point
        
        const_x = [1 findVal1 param.spline_length];
        
%     case 'startmid1endonknotabove'
%         mid1 = ceil(param.spline_length/2)+1;
%         [findVal1, findInd1] = findClosestValue(mid1, intermed_ind, 'above'); % find the closest knot point
%         
%         const_x = [1 findVal1 param.spline_length];        
        
    case 'startmid1endonknotbelow'
        mid1 = ceil(param.spline_length/2)+1;
        [findVal1, findInd1] = findClosestValue(mid1, intermed_ind, 'below'); % find the closest knot point
        
        const_x = [1 findVal1 param.spline_length];    

    case 'startmid2endonknot'
        mid1 = ceil(param.spline_length/3);
        [findVal1, findInd1] = findClosestValue(mid1, intermed_ind); % find the closest knot point
        mid2 = ceil(2*param.spline_length/3);
        [findVal2, findInd2] = findClosestValue(mid2, intermed_ind); % find the closest knot point
        
        const_x = [1 findVal1 findVal2 param.spline_length];
        
    case 'variableEdge'
        % find all the const_x in the window
        t_target = feature_full.t_cf_q_const;
        q_target = 1:length(t_target);
        const_x = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_target, q_target, param);
        const_x = [1 const_x param.spline_length];
        const_x = unique(const_x);
       
    case 'front'
end

switch cf_dq_constraint
    case 'constraint_dxddx'
        const_dx = const_x;
        const_ddx = const_x;
        
    case 'edge_dx'
        const_dx = [1 param.spline_length];
        const_ddx = [];
        
    case 'variable_dx'
        % find all the const_x in the window
        t_full = feature_full.t_cf_dq_const;
        q_full = feature_full.t_cf_dq_const;
        const_dx = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        
        const_ddx = [];
        
    case 'variableEdge_dx'
        % find all the const_x in the window
        t_full = feature_full.t_cf_dq_const;
        q_full = feature_full.t_cf_dq_const;
        const_dx = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        const_dx = [1 const_dx param.spline_length];
        
        const_ddx = [];        
        
    case 'variable_dxddx'
        % find all the const_x in the window
        t_full = feature_full.t_cf_dq_const;
        q_full = feature_full.t_cf_dq_const;
        const_dx = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        
        t_full = feature_full.t_cf_ddq_const;
        q_full = feature_full.t_cf_ddq_const;
        const_ddx = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);     
        
    case 'edge_dxddx'
        const_dx = [1 param.spline_length];
        const_ddx = [1 param.spline_length];
        
    case 'startmidend_dxddx'
        const_dx = [1 ceil(param.spline_length/2)+1 param.spline_length];
        const_ddx = [1 ceil(param.spline_length/2)+1 param.spline_length];
        
    case 'startmidmidend_dxddx'
        const_dx = [1 floor(param.spline_length/3) ceil(2*param.spline_length/3) param.spline_length];
        const_ddx = [1 floor(param.spline_length/3) ceil(2*param.spline_length/3) param.spline_length];
        
    case 'variableEdge_dxddx'
        % find all the const_x in the window
        t_full = feature_full.t_cf_dq_const;
        q_full = feature_full.t_cf_dq_const;
        const_dx = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        const_dx = [1 const_dx param.spline_length];
        const_dx = unique(const_dx);
        
        t_full = feature_full.t_cf_ddq_const;
        q_full = feature_full.t_cf_ddq_const;
        const_ddx = findMatch_t(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        const_ddx = [1 const_ddx param.spline_length];
        const_ddx = unique(const_ddx); 
        
    otherwise
        const_dx = [];
        const_ddx = [];
end

% insert the enforced constants so that we can calc the gradient properly
switch param.ioc_add_cf_as_knots
    case 'yes'
        intermed_ind = [intermed_ind const_x];
        
    case 'no'
        intermed_ind = [intermed_ind];
end

intermed_ind = unique(intermed_ind);

% add all the knot points as constraints
switch add_knot_as_cf_constraints
    case 'none'
        
    case 'yes'
        const_x = unique([intermed_ind const_x]); % this ensures there is a knot point at the end, and on constraints
end

%     const_x = const_x(1:3);

switch cf_q_constraint_y
    case 'variable'
        t_full = feature_full.t_cf_q_const;
        q_full = feature_full.y_cf_q_const;
        const_y = findMatch_y(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        
    case 'previous'
        const_y = feature_win.q(:, const_x); % the constraints are passed to the IOC and is compared against the normalized q, so should be normalized first
end

switch cf_dq_constraint_y
    case 'variable_dx'       
        t_full = feature_full.t_cf_dq_const;
        q_full = feature_full.y_cf_dq_const;
        const_dy = findMatch_y(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        const_ddy = [];
        
    case 'variable_dxddx'       
        t_full = feature_full.t_cf_dq_const;
        q_full = feature_full.y_cf_dq_const;
        const_dy = findMatch_y(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
        
        t_full = feature_full.t_cf_ddq_const;
        q_full = feature_full.y_cf_ddq_const;
        const_ddy = findMatch_y(t_win_findMatchArray, q_win_findMatchArray, t_full, q_full, param);
          
    case 'previous_dx'
        const_dy = feature_win.dq(:, const_dx);
        const_ddy = [];
        
    case 'previous_dxddx'
        const_dy = feature_win.dq(:, const_dx);
        const_ddy = feature_win.ddq(:, const_ddx);
        
    case 'zero_dx'
        const_dy = zero(size(feature_win.dq(:, const_dx)));
        const_ddy = [];
end

x_knot = param.t_spline(intermed_ind); % update the x knots to ensure the dt is correct
% x_knot = feature_win.t(intermed_ind);

param.intermed_ind_set = intermed_ind;
param.x_knot_set = x_knot;

param.const_x = const_x;
param.const_y = const_y;
param.const_dx = const_dx;
param.const_dy = const_dy;
param.const_ddx = const_ddx;
param.const_ddy = const_ddy;

% param.t_spline = feature_win.t;
param.dt_spline = mean(diff(feature_win.t));

spline_endcond_dx = [];
spline_endcond_dy = [];
spline_endcond_ddx = [];
spline_endcond_ddy = [];
spline_endcond_dddx = [];
spline_endcond_dddy = [];

% set the splining conditions so it matches the input data
switch sp_dq_endcond
    case 'constraint_dxddx'        
        for ind_const = 1:length(const_x)
            curr_constx = const_x(ind_const);
            spline_endcond_dx =  [spline_endcond_dx  curr_constx];
            spline_endcond_dy =  [spline_endcond_dy  feature_win.dq(:, curr_constx)];
            spline_endcond_ddx = [spline_endcond_ddx curr_constx];
            spline_endcond_ddy = [spline_endcond_ddy feature_win.ddq(:, curr_constx)];
%             spline_endcond_dddx = [spline_endcond_dddx curr_constx];
%             spline_endcond_dddy = [spline_endcond_dddy feature_win.dddq(:, curr_constx)];            
        end
        
    case 'knot_dxddx'
        for ind_const = 1:length(intermed_ind)
            curr_constx = intermed_ind(ind_const);
            spline_endcond_dx =  [spline_endcond_dx  curr_constx];
            spline_endcond_dy =  [spline_endcond_dy  feature_win.dq(:, curr_constx)];
            spline_endcond_ddx = [spline_endcond_ddx curr_constx];
            spline_endcond_ddy = [spline_endcond_ddy feature_win.ddq(:, curr_constx)];
%             spline_endcond_dddx = [spline_endcond_dddx curr_constx];
%             spline_endcond_dddy = [spline_endcond_dddy feature_win.dddq(:, curr_constx)];      
        end
        
        
    case 'edge_dx'
        %             if endcond_smoothedge_flag
%         dq_start =  feature_win.dq(:, 1)  ; % this data is passed directly to the spline and should not be normalized first
%         dq_end   =  feature_win.dq(:, end);
%         ddq_start = feature_win.ddq(:, 1) ;
%         ddq_end =   feature_win.ddq(:, end);

        spline_endcond_dx =  [1 param.spline_length];
        spline_endcond_dy =  feature_win.dq(:, [1 end]); 
%         spline_endcond_ddx = [1 param.spline_length];
%         spline_endcond_ddy = feature_win.ddq(:, [1 end]); 
%         spline_endcond_dddx = [1 param.spline_length];
%         spline_endcond_dddy = feature_win.dddq(:, [1 end]); 
        
        %               param.spline_endcond_y = zeros(param.dof_count, length(param.spline_endcond_x));
        
        case 'edge_dxddx'
        %             if endcond_smoothedge_flag
%         dq_start =  feature_win.dq(:, 1)  ; % this data is passed directly to the spline and should not be normalized first
%         dq_end   =  feature_win.dq(:, end);
%         ddq_start = feature_win.ddq(:, 1) ;
%         ddq_end =   feature_win.ddq(:, end);

        spline_endcond_dx =  [1 param.spline_length];
        spline_endcond_dy =  feature_win.dq(:, [1 end]); 
        spline_endcond_ddx = [1 param.spline_length];
        spline_endcond_ddy = feature_win.ddq(:, [1 end]); 
%         spline_endcond_dddx = [1 param.spline_length];
%         spline_endcond_dddy = feature_win.dddq(:, [1 end]); 
        
        %               param.spline_endcond_y = zeros(param.dof_count, length(param.spline_endcond_x));
        
    case 'variable_dxddx'
        % find all the const_x in the window
        const_x = [];
        for ind_const = 1:length(param.feature_full.t_sp_dq_const)
            currConstTime = param.feature_full.t_sp_dq_const(ind_const);
            [findVal, findInd] = findClosestValue(currConstTime, feature_win.t);
            if abs(findVal - currConstTime) < 0.99*param.dt_full % less than 1ms error
                const_x = [const_x findInd];
            end
        end
        
        for ind_const = 1:length(const_x)
            curr_constx = const_x(ind_const);
            spline_endcond_dx =  [spline_endcond_dx  curr_constx];
            spline_endcond_dy =  [spline_endcond_dy  feature_win.dq(:, curr_constx)];
            spline_endcond_ddx = [spline_endcond_ddx curr_constx];
            spline_endcond_ddy = [spline_endcond_ddy feature_win.ddq(:, curr_constx)];
%             spline_endcond_dddx = [spline_endcond_dddx curr_constx];
%             spline_endcond_dddy = [spline_endcond_dddy feature_win.dddq(:, curr_constx)];
        end
        
    case 'front_dxddx'
        curr_constx = const_x(1);
        spline_endcond_dx =  [spline_endcond_dx  curr_constx];
        spline_endcond_dy =  [spline_endcond_dy  feature_win.dq(:, curr_constx)];
        spline_endcond_ddx = [spline_endcond_ddx curr_constx];
        spline_endcond_ddy = [spline_endcond_ddy feature_win.ddq(:, curr_constx)];
%         spline_endcond_dddx = [spline_endcond_dddx curr_constx];
%         spline_endcond_dddy = [spline_endcond_dddy feature_win.dddq(:, curr_constx)];

    case 'none'
        
end

param.spline_endcond_dx = spline_endcond_dx;
param.spline_endcond_dy = spline_endcond_dy;
param.spline_endcond_ddx = spline_endcond_ddx;
param.spline_endcond_ddy = spline_endcond_ddy;
param.spline_endcond_dddx = spline_endcond_dddx;
param.spline_endcond_dddy = spline_endcond_dddy;


% a way to check the constraints
set_ioc_check_struct.intermed_ind_set = intermed_ind;
set_ioc_check_struct.x_knot_set = x_knot;

set_ioc_check_struct.const_x = const_x;
set_ioc_check_struct.const_y = const_y;
set_ioc_check_struct.const_dx = const_dx;
set_ioc_check_struct.const_dy = const_dy;
set_ioc_check_struct.const_ddx = const_ddx;
set_ioc_check_struct.const_ddy = const_ddy;

set_ioc_check_struct.dt_spline = mean(diff(feature_win.t));

% set_ioc_check_struct.spline_endcond_dx = spline_endcond_dx;
% set_ioc_check_struct.spline_endcond_dy = spline_endcond_dy;
% set_ioc_check_struct.spline_endcond_ddx = spline_endcond_ddx;
% set_ioc_check_struct.spline_endcond_ddy = spline_endcond_ddy;
% set_ioc_check_struct.spline_endcond_dddx = spline_endcond_dddx;
% set_ioc_check_struct.spline_endcond_dddy = spline_endcond_dddy;
set_ioc_check_struct.intermed_ind_setArray = set_ioc_check_struct.intermed_ind_set + startTime - 1;

set_ioc_check_struct
end

function q_out = findMatch_t(t_win, q_win, t_target, q_target, param)
% q_out = findMatch_t(feature_win.t, param.feature_full.t_cf_q_const, param.feature_full.t_cf_q_const, param)
    q_out = [];
    for ind_const = 1:length(t_target)
        currConstTime = t_target(ind_const);
        [findVal, findInd] = findClosestValue(currConstTime, t_win);
        if abs(findVal - currConstTime) < 0.99*param.dt_full % less than 1ms error
            q_out = [q_out q_win(findInd)];
        end
    end
    q_out = unique(q_out);

end

function q_out = findMatch_y(t_win, q_win, t_target, q_target, param)
% q_out = findMatch_t(feature_win.t, param.feature_full.t_cf_q_const, param.feature_full.t_cf_q_const, param)
    q_out = [];
    for ind_const = 1:length(t_win)
        currConstTime = t_win(ind_const);
        [findVal, findInd] = findClosestValue(currConstTime, t_target);
        if abs(findVal - currConstTime) < 0.99*param.dt_full % less than 1ms error
            q_out = [q_out q_target(:, findInd)];
        end
    end

end