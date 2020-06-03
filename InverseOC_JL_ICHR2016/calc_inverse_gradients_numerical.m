function [J_cost_all, J_const_all, cost_weight_coeff, param] = calc_inverse_gradients(feature_opt, c_gen, c_const, param)
% pull out knots from the traj, then generate a spline-based traj from that
y_opt_base{1} = feature_opt.q(:, param.intermed_ind);
y_opt_base{2} = feature_opt.dq(:, param.intermed_ind);
y_opt_base{3} = feature_opt.ddq(:, param.intermed_ind);

[y_opt_vec, param.matVec_struct] = mat2vec_multi(y_opt_base); % converting the state matrix (easier to read) to a vector (fmincon req)

len_ycombo = param.matVec_struct.feature_count*param.matVec_struct.feature_width*param.matVec_struct.entry_count;

% we don't want to use the actual c_gen array if it is passed in here. we
% just need it for its width
len_cost = size(c_gen, 2);
c_use = zeros(1, len_cost);

% J_h_array = cell(1, len_cost);
% g_h_array = cell(1, len_const);

% calculate the numerical values for the cost function
for ii_knot = 1:len_ycombo
for ii_dofs = 1
for ii_h = 1:length(param.h_array)
    % iterate through all the states (knot/dof combination) and generate
    % the values we need for the numerical gradient calculations
    
    % calculate the spline corresponding to this knot/dof/h combination
    curr_spline_fit.q = applyKnotOffsetAndSpline(y_opt_vec, feature_opt, ii_knot, ii_h, param);

    % using the spline, calculate the features
    feature_use = calc_features(curr_spline_fit, [], param);
    
    % compute the cost function values corresponding to this combination.
    % that is, how does the cost function change when we change the
    [~, ~, J_cost_array] = calc_direct_cost(c_use, feature_use, param);
    
    % also the constraint functions corresponding
    [ceq, ceq_x, ceq_dx, ceq_ddx] = calc_direct_const(c_const, feature_use, param);
    len_const = length(ceq);
    
    % now save the variable data that we need it later. normalize the data
    % before inserting it into the h array
%     feature_use = normalize_features(feature_use, param);
   
%     % the dim of this 5D array: DOF, time, h offset, knot that was shifted,
%     % time entry that was shifted. note that these values are currently not
%     % normalized as it has been commented out from above
%     param.q_h_array(:, :, ii_h, ii_knot, ii_dofs) =     feature_use.q;
%     param.dq_h_array(:, :, ii_h, ii_knot, ii_dofs) =    feature_use.dq;
%     param.ddq_h_array(:, :, ii_h, ii_knot, ii_dofs) =   feature_use.ddq;
%     param.dddq_h_array(:, :, ii_h, ii_knot, ii_dofs) =  feature_use.dddq;
%     param.x_h_array(:, :, ii_h, ii_knot, ii_dofs) =     feature_use.x;
%     param.dx_h_array(:, :, ii_h, ii_knot, ii_dofs) =    feature_use.dx;
%     param.ddx_h_array(:, :, ii_h, ii_knot, ii_dofs) =   feature_use.ddx;
%     param.dddx_h_array(:, :, ii_h, ii_knot, ii_dofs) =  feature_use.dddx;
%     param.tau_h_array(:, :, ii_h, ii_knot, ii_dofs) =   feature_use.tau;
%     param.dtau_h_array(:, :, ii_h, ii_knot, ii_dofs) =  feature_use.dtau;
%     param.ddtau_h_array(:, :, ii_h, ii_knot, ii_dofs) = feature_use.ddtau;
%     param.cop_h_array(:, :, ii_h, ii_knot, ii_dofs) =   feature_use.cop;
%     param.dcop_h_array(:, :, ii_h, ii_knot, ii_dofs) =  feature_use.dcop;
%     param.ddcop_h_array(:, :, ii_h, ii_knot, ii_dofs) = feature_use.ddcop;
%     param.com_h_array(:, :, ii_h, ii_knot, ii_dofs) =   feature_use.com;
%     param.dcom_h_array(:, :, ii_h, ii_knot, ii_dofs) =  feature_use.dcom;
%     param.ddcom_h_array(:, :, ii_h, ii_knot, ii_dofs) = feature_use.ddcom;
%     param.ep_h_array(:, :, ii_h, ii_knot, ii_dofs) =    feature_use.ep;
%     param.ek_h_array(:, :, ii_h, ii_knot, ii_dofs) =    feature_use.ek;
%     param.geo_h_array(:, :, ii_h, ii_knot, ii_dofs) =   feature_use.geo;
%     param.en_h_array(:, :, ii_h, ii_knot, ii_dofs) =    feature_use.en;
    
    % the order of the spline will be a bit different. the dim of these
    % entries are dof, knot, h, for quicker access than the above entries
    param.spline_fit(ii_dofs, ii_knot, ii_h) = curr_spline_fit;
    
    % same idea goes for the cost function, but we can make the cost
    % function more automatic instead of having to modify it all the time
    % by loading them into cell arrays
    for ii_J = 1:len_cost
        J_h_array{ii_J}(ii_dofs, ii_knot, ii_h) = J_cost_array(ii_J);
    end
    
    % and now build the constraint functions
    for ii_J = 1:len_const
        g_h_array{ii_J}(ii_dofs, ii_knot, ii_h) = ceq(ii_J);
    end
    
    if 0 % if we want to check the quality of the gradient
        y_opt = feature_opt.q';
        dy = calcDeriv(y_opt, param.dt)';

        feature_unnorm = unnormalize_features(feature_use, param);
        q_unnorm = feature_unnorm.q;
        dq_unnorm = feature_unnorm.dq;

        figure; 
        subplot(211); plot(y, 'x'); hold on; plot(q_unnorm')
        plot(param.intermed_ind, y(param.intermed_ind, :), 'ko'); 

        subplot(212); plot(dy, 'x'); hold on; plot(dq_unnorm')
        plot(param.intermed_ind, dy(param.intermed_ind, :), 'ko');
%             figure; plot(ddy, 'x'); hold on; plot(ddq')
    end
end
end
end

% % % % plot the generated plot vs the spline
%     plot_inverse;

% use the calculated values to generate the proper cost functions and
% constraint matrix. the length order is:
%   t1, x1
%   t1, x2
%   t1, x3
%   t2, x1
%   t2, x2
%   t2, x3 ... etc
J_template = zeros(len_ycombo, 1);

gradJ_cost = cell(1, len_cost);
gradJ_cost(:) = {J_template}; % preallocate the cost function array

gradJ_const = cell(1, len_const);
gradJ_const(:) = {J_template}; % preallocate the cost function array

for ii_knot = 1:len_ycombo
for ii_dofs = 1
    
%     ind_J = (ii_knot-1)*param.dof_count + ii_dofs;
    ind_J = (ii_knot);
    
    % iterate through all th entries in the J array, as calculated above,
    % to note how they change corresponding to the changes in q_knots
    for ii_J = 1:len_cost
        J_curr_knotdofs = J_h_array{ii_J}(ii_dofs, ii_knot, :); % this value is normalized since the normalization function is in the cost calculations
        gradJ_cost{ii_J}(ind_J) = calc_diff_num(J_curr_knotdofs, param, param.diffType, 'cost');
    end
    
    % need to break the full length array into 3 smaller components to
    % check which (ie x, dx or ddx) is currently active to determine where
    % the constraints are
    for ii_J = 1:len_const
        J_curr_knotdofs = g_h_array{ii_J}(ii_dofs, ii_knot, :); % this value is normalized since the normalization function is in the cost calculations
        gradJ_const{ii_J}(ind_J) = calc_diff_num(J_curr_knotdofs, param, param.diffType, 'cost');
    end
end
end

% combine these values into the least squares form, for combining later on
% J_cost_all = {gradJ_cost_ddq, gradJ_cost_ddx, gradJ_cost_tau};
% [gradJ_cost{:}] 
% [gradJ_const{:}]
J_cost_all = gradJ_cost;
J_const_all = gradJ_const;

% pull out the actual cost function. at h=0, it shouldn't have shifted 
cost_weight_coeff = zeros(1, len_cost);
h_ind1 = param.h_vals == 0;
for ind = 1:len_cost
    cost_weight_coeff(ind) = J_h_array{ind}(1, 1, h_ind1); 
end
end

function gradJ_const_x = matchKnotsForConstraints(ii_knot, ii_dofs, ind_J, const_x, const_y, q_knotdofs, gradJ_const_x, param)
    knotVal = param.intermed_ind(ii_knot);
    const_ind = find(knotVal == const_x);

    if ~isempty(const_ind)
        extraParam.knotVal = knotVal;
        extraParam.dofVal = ii_dofs;
        extraParam.y_base = const_y(ii_dofs, const_ind);
        diffVal = calc_diff_num(q_knotdofs, param, param.diffType, 'const', extraParam);
        
        if diffVal == 0
            gradJ_const_x(ind_J) = 1e-30; % it's a constraint that did not change wrt to q
        else
            gradJ_const_x(ind_J) = diffVal;
        end
        
    end
end