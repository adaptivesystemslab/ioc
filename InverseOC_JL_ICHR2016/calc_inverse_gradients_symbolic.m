function [J_cost_all, J_const_all, param] = calc_inverse_gradients(feature_opt, c_gen, c_const, param)
% pull out knots from the traj, then generate a spline-based traj from that
y_opt_base{1} = feature_opt.q(:, param.intermed_ind);
y_opt_base{2} = feature_opt.dq(:, param.intermed_ind);
y_opt_base{3} = feature_opt.ddq(:, param.intermed_ind);
curr_spline_fit.q  =  calc_direct_spline(y_opt_base, param);
curr_spline_fit.syms = 0;
feature_use_val = calc_features(curr_spline_fit, feature_opt, param);

curr_spline_fit.syms = 1;
feature_use_syms = calc_features_syms(curr_spline_fit, [], param);
y_opt_syms{1} = feature_use_syms.q(:, param.intermed_ind);
y_opt_syms{2} = feature_use_syms.dq(:, param.intermed_ind);
y_opt_syms{3} = feature_use_syms.ddq(:, param.intermed_ind);
[y_opt_vec, param.matVec_struct] = mat2vec_multi_syms(y_opt_syms); % converting the state matrix (easier to read) to a vector (fmincon req)

 
% compute the cost function values corresponding to this combination.
% that is, how does the cost function change when we change the
cost_len = size(c_gen, 2);
c_use = zeros(1, cost_len);
[~, ~, J_h_array] = calc_direct_cost_syms(c_use, feature_use_syms, param);

x_knot_ind = 1:length(param.x_knot);
len_x_knot_ind = length(x_knot_ind);
J_template = zeros(param.dof_count*len_x_knot_ind, 1);

% gradJ_cost = cell(1, cost_len);
% gradJ_cost(:) = {J_template}; % preallocate the cost function array
% gradJ_const = cell(1, 3);
% gradJ_const(:) = {J_template}; % preallocate the const function array

% now assemble the constraint function
gradJ_const_x = sym(zeros(length(param.x_knot)*param.dof_count, 1));
gradJ_const_dx = sym(zeros(length(param.x_knot)*param.dof_count, 1));
gradJ_const_ddx = sym(zeros(length(param.x_knot)*param.dof_count, 1));

len_cost = size(c_gen, 2);

len_ycombo = param.matVec_struct.feature_count*param.matVec_struct.feature_width*param.matVec_struct.entry_count;

h_h_array = calc_direct_const_sym(c_const, feature_use_syms, param);
len_const = size(h_h_array, 1);

for ii_J = 1:len_const
    for ii_knot = 1:len_ycombo
        % assemble the function
        J_curr_knotdofs = diff(h_h_array(ii_J), y_opt_vec(ii_knot)); % this value is normalized since the normalization function is in the cost calculations
        gradJ_const{ii_J}(ii_knot) = J_curr_knotdofs;
    end
end

for ii_knot = 1:len_ycombo
    for ii_dofs = 1:param.dof_count
        % load the appropriate slice of variables
        ind_J = (ii_knot-1)*param.dof_count + ii_dofs;
        
        % diff the cost functions
        for ii_J = 1:len_cost
            J_curr_knotdofs = diff(J_h_array{ii_J}, y_opt_vec(ii_knot)); % this value is normalized since the normalization function is in the cost calculations
            gradJ_cost{ii_J}(ind_J) = J_curr_knotdofs;
        end
    end
end

% need to figure out which stage needs diffing and which stage can be saved







param.h_h_array{1} = gradJ_const_x; 
param.h_h_array{2} = gradJ_const_dx; 
param.h_h_array{3} = gradJ_const_ddx; 

% pull out the denominators for the gradient calculation
gradSyms_q = mat2vec(feature_use_syms.q(:, param.intermed_ind(:)));

% iterate through all th entries in the J array, as calculated above,
% to note how they change corresponding to the changes in q_knots
for ii_J = 1:length(param.J_h_array)
    J_curr_knotdofs = param.J_h_array{ii_J}; % this value is normalized since the normalization function is in the cost calculations
    J_curr_diffed = gradient(J_curr_knotdofs, gradSyms_q); 
    J_curr_subbed = sym2val(J_curr_diffed, feature_use_syms, feature_use_val);
    
    gradJ_cost{ii_J} = J_curr_subbed;
end

for ii_J = 1:length(param.h_h_array)
    J_curr_knotdofs = param.h_h_array{ii_J}; % this value is normalized since the normalization function is in the cost calculations
    J_curr_diffed = jacobian(J_curr_knotdofs, gradSyms_q); 
    J_curr_diag = diag(J_curr_diffed);
    J_curr_subbed = sym2val(J_curr_diag, feature_use_syms, feature_use_val);
    
    gradJ_const{ii_J} = J_curr_subbed;
end

% combine these values into the least squares form, for combining later on
% J_cost_all = {gradJ_cost_ddq, gradJ_cost_ddx, gradJ_cost_tau};
J_cost_all = gradJ_cost;
J_const_all = gradJ_const;

% [gradJ_cost{:}]
% [gradJ_const{:}]
end

function constFct = matchKnotsForConstraints(ii_knot, ii_dofs, const_x, const_y, varBase_syms, param)
    knotVal = param.intermed_ind(ii_knot);
    const_ind = find(knotVal == const_x);

    if ~isempty(const_ind)
        constFct = varBase_syms(ii_dofs, knotVal) - const_y(ii_dofs, const_ind);
    else
        constFct = 0;
    end
end

function diffVal = matchKnotsForConstraints_chainRule(ii_knot, ii_dofs, const_x, const_y, varBase_syms, varDiff_syms, dq_syms, param)
%     % calculate the specific entries for the constraint
%     if isfield(param, 'const_x')
%         gradJ_const_x(ind_J) = matchKnotsForConstraints_chainRule(ii_knot, ii_dofs, param.const_x, param.const_y, ...
%             feature_use_syms.q, feature_use_syms.dq, feature_use_syms.dq, param);
%         
%     end
%     
%     if isfield(param, 'const_dx')
%         gradJ_const_dx(ind_J) = matchKnotsForConstraints_chainRule(ii_knot, ii_dofs, param.const_dx, param.const_dy, ...
%             feature_use_syms.dq, feature_use_syms.ddq, feature_use_syms.dq, param);
%     end
%     
%     if isfield(param, 'const_ddx')
%         gradJ_const_ddx(ind_J) = matchKnotsForConstraints_chainRule(ii_knot, ii_dofs, param.const_ddx, param.const_ddy, ...
%             feature_use_syms.ddq, feature_use_syms.dddq, feature_use_syms.dq, param);
%     end

    knotVal = param.intermed_ind(ii_knot);
    const_ind = find(knotVal == const_x);

    if ~isempty(const_ind)
        % the following is the constraint function. since we are going to
        % diff it, the constant will drop away, and we will work with just
        % the starting part
%         diffVal = q_syms(ii_dofs, knotVal) - const_y(ii_dofs, const_ind);
        
        % J = dq - dq_set
        % dJ/dq1 = d(dq1)/d(q1) = d(dq1)/d(t) * d(t)/d(q1) = ddq1 * dq1^-1
        base_diff = varDiff_syms(ii_dofs, knotVal);
        dq_diff = dq_syms(ii_dofs, knotVal);
        diffVal = base_diff * (1/(dq_diff));        
    else
        diffVal = sym(0);
    end
end