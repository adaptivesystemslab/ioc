function [c_out_pivot, J_out_contrib_ratio, output_inverse, cost_functions_used, J_cost_use_combined, J_costconst_check] = ...
    main_inverse(feature_opt, c_gen, c_const, param, param_opt, c_init_guess)

% 3dof model squat simulation
% inverse optimization - given trajectory, we want to recover cost function

displaySetting = 'off'; % 'final' 'notify' 'off' 'iter-detailed'

param = update_intermed_ind(param, param_opt.knot_count);

param.h_vals =  [0 -1/2 +1/2];
param.h_array = param.h_vals*param.h; 

switch param.ioc_gradient_method
    case 'numerical'
        [J_cost_num, J_const_num, param] = calc_inverse_gradients_numerical(feature_opt, c_gen, c_const, param);
        J_cost_all = J_cost_num;
        J_const_all = J_const_num;
        
    case 'analytical_spline'
        [J_cost_all, J_const_all, param] = calc_inverse_gradients_analytical(feature_opt, c_gen, param);
        
    case 'symbolic'
        [J_cost_sym, J_const_sym, param] = calc_inverse_gradients_symbolic(feature_opt, c_gen, param);
        J_cost_all = J_cost_sym;
        J_const_all = J_const_sym;
end

J_costconst_check = [J_cost_num{:} J_const_num{:}];

if 0    
    cost_const_num = [J_cost_num{:} J_const_num{:}]
    cost_const_sym = [J_cost_sym{:} J_const_sym{:}]
end

J_cost_check = J_cost_all;

J_cost_use_combined = horzcat(J_cost_check{:}); 
corrValArray = corrcoef(J_cost_use_combined);
corrValArray_full = corrValArray;

[~                  , corrValArray_dotp111, corrValArrayOut_dotp111]  = checkCorrelation(J_cost_all, J_const_all, 'dotp_111', param);
[~                  , corrValArray_dotp100, corrValArrayOut_dotp111]  = checkCorrelation(J_cost_all, J_const_all, 'dotp_100', param);
[~                  , corrValArray_corr_111, corrValArrayOut_corr111] = checkCorrelation(J_cost_all, J_const_all, 'corr_111', param);
[~                  , corrValArray_corr_100, corrValArrayOut_corr100] = checkCorrelation(J_cost_all, J_const_all, 'corr_100', param);

corrJHqRangeMtx = corrValArray_dotp100(1:3, 4:end);
corrJHqmax = max(max(corrJHqRangeMtx));
corrJHqmin = min(min(corrJHqRangeMtx));
corrJHqrange = corrJHqmax - corrJHqmin;

[cost_functions_used, corrValArray, corrValArrayOut, corrValId] = checkCorrelation(J_cost_all, {}, 'corr_111', param);

cost_functions_used = cost_functions_used(cost_functions_used > 0);
J_cost_use = J_cost_check(cost_functions_used);
param.cost_functions_ioc = cost_functions_used; % save this for use inside the cost function calculation
len_cost_functions_full = length(c_gen);
len_cost_functions_used = length(cost_functions_used);
c_init_used = c_init_guess(cost_functions_used);

% create the pivot array. this is used to iterate through all the potential
% pivot values
fullPivotSet = cell(1, len_cost_functions_used - 1);
for ind_pivotCount = 1:len_cost_functions_used - 1
    fullPivotSet{ind_pivotCount} = ['c_use(' num2str(ind_pivotCount) ')'];
end

ind_pivotToUse = 0;
for ind_pivotset = 1:len_cost_functions_full
    % the unknown parameters are...
    % lambda: 3 dof x 3 each (q_init, q_mid, and q_final) = 21 constraints
    % c: 2 cost function coeffs            
    if isempty(find(cost_functions_used == ind_pivotset, 1))
        % the entry we're on right now is not part of the unique set. we
        % don't want to run the IOC but need to create the rest of the
        % variables to keep the indices consistent
        
        lambda_init = zeros(param.dof_count, length(param.const_x) + length(param.const_dx) + length(param.const_ddx));
        c_init = 1*ones(len_cost_functions_used - 1, 1);
        [z_init, param.feature_count, param.c_array_len] = mat2vec_extra(lambda_init, c_init);
        
        J_cost_all_merge = horzcat(J_cost_all{:});
        size_J_cost_merge = size(J_cost_all_merge);
        J_coeff_out = zeros(size_J_cost_merge); 
        J_coeff_total_out = zeros(size(J_cost_all{1})); 
        resnorm = 1e6;
        EXITFLAG = -1;
        A = zeros(length(z_init), length(z_init));
        gradDotSign = [];
        resnorm_lsqlin = 0;
        residual_lsqlin = 0;
        J_out_contrib{ind_pivotset} = zeros(1, len_cost_functions_full);
        J_out_contrib_ratio{ind_pivotset} = zeros(1, len_cost_functions_full);
        J_out_contrib_array{ind_pivotset} = zeros(1, len_cost_functions_full);
        c_out_pivot{ind_pivotset} = zeros(1, len_cost_functions_full);
        J_out_contrib_original{ind_pivotset} = zeros(1, len_cost_functions_full);
        c_out_pivot_original{ind_pivotset} = zeros(1, len_cost_functions_full);
        lambda_out{ind_pivotset} = mat2vec(lambda_init);
        
    else
        % otherwise, run the set    
        ind_pivotToUse = ind_pivotToUse + 1;

        param.pivot = [fullPivotSet(1:ind_pivotToUse-1) '1' fullPivotSet(ind_pivotToUse:end)];
        param.pivotTerm = ind_pivotToUse;
        
        iterLoop = 1;
        lambda_set = ones(1, length(J_const_all));
        c_set = [c_init_used(1:ind_pivotToUse-1) c_init_used(ind_pivotToUse+1:end)]'; % use the initial guess passed in from the outside
                    Mesoptions = optimset('TolFun',param_opt.TolFun,'TolX',param_opt.TolX,'TolCon',param_opt.TolCon,...
                'display',displaySetting,'MaxIter',param_opt.maxiter,'LargeScale','off','Diagnostics','off');
                
        for loopInd = 1:10
            % run for lambda           
            [z_init, param.matVec_struct] = mat2vec_multi([], c_set);

            [~, linConst_C, ~, lsqlin_d] = ioc_cost_func3(z_init, lambda_set, J_cost_use, J_const_all, c_const, param);
            lsqlin_d = -lsqlin_d; % this is to get Ax+b, since we're doing Ax-(-b)
            lsqlin_A = zeros(size(linConst_C, 2), size(linConst_C, 2));
            lsqlin_A(1:length(c_set), 1:length(c_set)) = diag(-ones(size(c_set))); % the constrains are Ax<b, so A should be neg for the cost weights
            lsqlin_b = zeros(size(z_init));            
            lsqlin_Aeq = [];
            lsqlin_beq = [];
            
            [z_calc,resnorm_lsqlin,residual_lsqlin,EXITFLAG,output_lsqlin] = lsqlin(linConst_C,lsqlin_d,lsqlin_A,lsqlin_b,lsqlin_Aeq,lsqlin_beq,[],[],z_init,Mesoptions);
            c_set = z_calc;
            
            % run for c
            [z_init, param.matVec_struct] = mat2vec_multi(lambda_set, []);
            
            [~, linConst_C, ~, lsqlin_d] = ioc_cost_func4(z_init, c_set, J_cost_use, J_const_all, c_const, param);
            lsqlin_d = -lsqlin_d; % this is to get Ax+b, since we're doing Ax-(-b)
            lsqlin_A = zeros(size(linConst_C, 2), size(linConst_C, 2));
            lsqlin_A(1:length(c_set), 1:length(c_set)) = diag(-ones(size(c_set))); % the constrains are Ax<b, so A should be neg for the cost weights
            lsqlin_b = zeros(size(z_init));            
            lsqlin_Aeq = [];
            lsqlin_beq = [];
            
            [z_calc,resnorm_lsqlin,residual_lsqlin,EXITFLAG,output_lsqlin] = lsqlin(linConst_C,lsqlin_d,lsqlin_A,lsqlin_b,lsqlin_Aeq,lsqlin_beq,[],[],z_init,Mesoptions);            
        end
        
        % revert the state vector back into array form so we can extract the
        % recovered weights
        [lambda_cell, c_use] = vec2mat_multi(z_calc, param.matVec_struct);
        lambda = lambda_cell{1};
        lambda_out{ind_pivotset} = mat2vec(lambda);
        c_out_temp = zeros(size(param.pivot));
        for ind_innerset = 1:length(param.pivot) % copy out the weights
            c_out_temp(ind_innerset) = eval(param.pivot{ind_innerset});
        end

        c_out_pivot_local = zeros(1, len_cost_functions_full);
        c_out_pivot_local(cost_functions_used) = c_out_temp; % calculate the original 
        c_out_pivot_original{ind_pivotset} = c_out_pivot_local; % re-align the cost functions by zero'ing the dropped colinear entries
        
        % calculate the contribution ratio of the recovered weights
        h_ind1 = param.h_vals == 0;
        y_opt_base = feature_opt.q(:, param.intermed_ind);
        splineFit.q = calc_direct_spline(y_opt_base, param);
        feature_use = calc_features(splineFit, [], param);

        [~, J_out_contrib_original{ind_pivotset}] = calc_direct_cost(c_out_pivot_original{ind_pivotset}, feature_use, param);
%         [J, A, x, b, J_coeff_out] = ioc_cost_func(z_calc, J_cost_use, J_const_all, c_const, param);

        % now apply the zeroing and remove entries that are below the
        % threshold
        c_out_pivot_zeroing = c_out_pivot_local;
        J_out_pivot_zeroing = J_out_contrib_original{ind_pivotset};
        c_out_pivot_zeroing(abs(c_out_pivot_zeroing) < param.recoverThreshold_c) = 0;
        c_out_pivot_zeroing(abs(J_out_pivot_zeroing) < param.recoverThreshold_j) = 0;
        c_out_pivot{ind_pivotset} = c_out_pivot_zeroing; % re-align the cost functions by zero'ing the dropped colinear entries
        [J_out_contrib{ind_pivotset}, J_out_contrib_ratio{ind_pivotset}, J_out_contrib_array{ind_pivotset}] = ...
            calc_direct_cost(c_out_pivot{ind_pivotset}, feature_use, param);
        
        c_zeroed = [c_out_pivot_zeroing(1:ind_pivotToUse-1) c_out_pivot_zeroing(ind_pivotToUse+1:end)]'; % use the initial guess passed in from the outside
        [z_zeroed, param.matVec_struct] = mat2vec_multi(lambda, c_zeroed);
        [J, A, x, b, J_coeff_out] = ioc_cost_func2(z_zeroed, J_cost_use, J_const_all, c_const, param);
        
        resnorm = norm(A*x+b)/length(b);

        % values to look at
        J_check = J_costconst_check;
        J_out_contrib_array_look = J_out_contrib_array{ind_pivotset};
        
        weight_look = c_out_pivot{ind_pivotset};
        z_calclook = z_calc(3:end)';
        
        Ax = (A*x)';
        b_look = b';
        AXB_look = [Ax; b_look];
        J_ioc_look = [J' sum(J) norm(J)/length(b)];
        
         J_out_contriub_look = J_out_contrib{ind_pivotset};
%         ind_pivotset = 4; c_out_pivot{ind_pivotset} = [1 1 1];
%           resnorm_look = norm(A*x+b)/length(b);
        
        % calculate the sign of the dot prod between the pivot and all others
        gradDotSign = [];
        for ind_cgen1 = 1:len_cost_functions_used
    %         for ind_cgen2 = 1:length(c_gen)
                dotProd = dot(b, J_coeff_out(:, ind_cgen1));
                gradDotSign(ind_cgen1) = sign(dotProd);
    %         end
        end

        gradDotSign(ind_pivotset) = 0; % ignore dot prod of the pivot against itself

        c_gen_mod = c_gen; c_gen_mod(c_gen == 0) = 1; % replace the version we're saving so we don't multiply the cost weights by 0

        J_coeff_total_out = 0;
        for ind_pivot = 1:length(c_out_temp)
            J_coeff_total_out = J_coeff_total_out + c_out_temp(ind_pivot)*J_coeff_out(:, ind_pivot);
        end

        % if we removed entries from the pivot before, reshape the matrix so
        % it's properly shaped
    %     J_out_contrib_local = zeros(1, len_cost_functions_full);
    %     J_out_contrib_local(cost_functions_used) = J_out;
    end
    
    % saving the gradient output
    output_inverse(ind_pivotset).J_coeff_individual_out = J_coeff_out;
    output_inverse(ind_pivotset).J_coeff_total_out = J_coeff_total_out;
    output_inverse(ind_pivotset).resnorm = resnorm; 
    output_inverse(ind_pivotset).exitflag = EXITFLAG;
    output_inverse(ind_pivotset).A = A;
    output_inverse(ind_pivotset).condA = cond(A);
    output_inverse(ind_pivotset).svd_s = svd(A);
    output_inverse(ind_pivotset).gradDotSign = gradDotSign;
    output_inverse(ind_pivotset).corrValMatrix_All = corrValArray_full;
    output_inverse(ind_pivotset).corrValMatrix_NaN = corrValArray;
    output_inverse(ind_pivotset).corrValArray = corrValArrayOut';
    output_inverse(ind_pivotset).corrValId = corrValId;
    output_inverse(ind_pivotset).J_out_array = J_out_contrib_array{ind_pivotset};
    output_inverse(ind_pivotset).J_out_contrib = J_out_contrib_ratio{ind_pivotset};
    output_inverse(ind_pivotset).c_recovered = c_out_pivot{ind_pivotset};
    output_inverse(ind_pivotset).c_offset = c_out_pivot{ind_pivotset}*c_gen_mod(ind_pivotset);
    output_inverse(ind_pivotset).rmse = 0; % this will be updated outside the fct
    output_inverse(ind_pivotset).c_out_pivot_original = c_out_pivot_original{ind_pivotset}; 
    output_inverse(ind_pivotset).J_out_contrib_original = J_out_contrib_original{ind_pivotset};
    output_inverse(ind_pivotset).lambda = lambda_out{ind_pivotset};
    output_inverse(ind_pivotset).resnorm_lsqlin = resnorm_lsqlin;
    output_inverse(ind_pivotset).residual_lsqlin = residual_lsqlin;
    output_inverse(ind_pivotset).corrJHqmax = corrJHqmax;
    output_inverse(ind_pivotset).corrJHqmin = corrJHqmin;
    output_inverse(ind_pivotset).corrJHqrange = corrJHqrange;
end

% if exist('S', 'var')
%     format long
%     c_list = c_gen;
%     for ind_innerset = 1:length(param.pivot)
%         c_list = [c_list; c_out_pivot{ind_innerset}*c_gen(ind_innerset)];
%     end
%     
%     c_list
% end
end