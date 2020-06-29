function [c_out_pivot_normalized, J_out_contrib_ratio, output_inverse, cost_functions_used_correlation, J_cost_use_combined, J_costconst_check] = ...
    main_inverse(feature_opt, c_gen, c_const, param, param_opt, c_init_guess, previousUniqueFunctions, saveFolder)

% 3dof model squat simulation
% inverse optimization - given trajectory, we want to recover cost function

displaySetting = 'off'; % 'final' 'notify' 'off' 'iter-detailed'

param = update_intermed_ind(param, param_opt.knot_count);

param.h_vals =  [0 -1/2 +1/2];
param.h_array = param.h_vals*param.h; 

switch param.ioc_gradient_method
    case 'numerical'
        [J_cost_num, J_const_num, cost_weight_coeff, param] = calc_inverse_gradients_numerical(feature_opt, c_gen, c_const, param);
        J_cost_all = J_cost_num;
        J_const_all = J_const_num;
        
%     case 'analytical_spline'
%         [J_cost_all, J_const_all, param] = calc_inverse_gradients_analytical(feature_opt, c_gen, param);
        
    case 'symbolic'
        [J_cost_sym, J_const_sym, param] = calc_inverse_gradients_symbolic(feature_opt, c_gen, c_const, param);
        J_cost_all = J_cost_sym;
        J_const_all = J_const_sym;
end

J_costconst_check = [J_cost_num{:} J_const_num{:}];

% if 0    
%     cost_const_num = [J_cost_num{:} J_const_num{:}]
%     cost_const_sym = [J_cost_sym{:} J_const_sym{:}]
% end

J_cost_check = J_cost_all;
J_zero = zeros(size(J_cost_all{1}));
J_cost_use_combined = horzcat(J_cost_check{:}); 
cost_functions_used_full = 1:length(J_cost_check);
J_cost_use = J_cost_check;
len_cost_functions_full = length(c_gen);

y_opt_base{1} = feature_opt.q(:, param.intermed_ind);
y_opt_base{2} = feature_opt.dq(:, param.intermed_ind);
y_opt_base{3} = feature_opt.ddq(:, param.intermed_ind);
splineFit.q = calc_direct_spline(y_opt_base, param);
feature_use_costfctCalc = calc_features(splineFit, [], param);

if param.removeHighlyCorrelatedJ
    [J_out_contrib_full, J_out_contrib_ratio_full, J_out_contrib_array_full] = calc_direct_cost(ones(1, len_cost_functions_full), feature_use_costfctCalc, param);
    [cost_functions_used_correlation, rref_out_full, eig_vals_full, eig_cutoff_full] = checkCorrelationSVD(J_cost_use_combined, param, J_out_contrib_array_full, 'prevWin');
else
    cost_functions_used_correlation = cost_functions_used_full;
    rref_out_full = [];
    eig_vals_full = [];
    eig_cutoff_full = [];
end

J_cost_use_reduced = {};
% zero out the pivots that are correlated
for ind_pivotset = 1:len_cost_functions_full
    if isempty(find(cost_functions_used_correlation == ind_pivotset, 1))
        J_cost_use{ind_pivotset} = J_zero; 
    else
        J_cost_use_reduced{end+1} = J_cost_use{ind_pivotset};
    end
end

param.cost_functions_ioc = cost_functions_used_full; % save this for use inside the cost function calculation
len_cost_functions_used_correlation = length(cost_functions_used_correlation);
cost_functions_used_correlation_array = 1:length(cost_functions_used_correlation);
cost_function_ind = 0;
c_init_used = c_init_guess(cost_functions_used_correlation);

% create the pivot array. this is used to iterate through all the potential
% pivot values
fullPivotSet = cell(1, len_cost_functions_used_correlation - 1);
for ind_pivotCount = 1:len_cost_functions_used_correlation - 1
    fullPivotSet{ind_pivotCount} = ['c_use(' num2str(ind_pivotCount) ')'];
end

switch param.pivotSelection
    case 'all'
        % don't remove any cost functions
        cost_function_used_pivot = cost_functions_used_correlation;
        
%     case 'minJ'
%         cost_weight_coeff(
%         [minVal, minInd] = min(cost_weight_coeff);
%         cost_function_used_pivot = minInd;
end

% ind_pivotToUse = 0;
for ind_pivotset = 1:len_cost_functions_full
    % the unknown parameters are...
    % lambda: 3 dof x 3 each (q_init, q_mid, and q_final) = 21 constraints
    % c: 2 cost function coeffs            
% % %     if isempty(find(cost_functions_used_correlation == ind_pivotset, 1))
    if isempty(find(cost_function_used_pivot == ind_pivotset, 1))
        % the entry we're on right now is not part of the unique set. we
        % don't want to run the IOC but need to create the rest of the
        % variables to keep the indices consistent
        
        lambda_init = zeros(1, length(J_const_num));
        c_init = 1*ones(len_cost_functions_used_correlation - 1, 1);
        [z_init, param.feature_count, param.c_array_len] = mat2vec_extra(lambda_init, c_init);
        
        J_cost_all_merge = horzcat(J_cost_all{:});
        size_J_cost_merge = size(J_cost_all_merge);
        J_coeff_out = zeros(size_J_cost_merge); 
        J_coeff_total_out = zeros(size(J_cost_all{1})); 
        resnorm_constrained = 1e6;
        EXITFLAG = -1;
        A = zeros(size(J_cost_all_merge, 1), size(J_cost_all_merge, 1));
        gradDotSign = [];
        c_normFactor = 1;
        resnorm_lsqlin_orig = 1e6;
        resnorm_lsqlin = 1e6;
        residual_lsqlin = 1e6;
        resnorm_lsqlin_unconst = 1e6;
        c_gen_mod = zeros(size(c_gen));
        c_out_temp_unconst = zeros(size(c_gen));
        J_out_contrib{ind_pivotset} = 1e6;
        J_out_contrib_ratio{ind_pivotset} = zeros(1, len_cost_functions_full);
        J_out_contrib_array{ind_pivotset} = zeros(1, len_cost_functions_full);
        c_out_pivot{ind_pivotset} = zeros(1, len_cost_functions_full);
%         J_out_contrib_original{ind_pivotset} = zeros(1, len_cost_functions_used_correlation);
        c_out_pivot_original{ind_pivotset} = zeros(1, len_cost_functions_full);
        c_out_pivot_normalized{ind_pivotset} = zeros(1, len_cost_functions_full);
        lambda_out{ind_pivotset} = mat2vec(lambda_init);        
    else
        % otherwise, run the set 
        cost_function_ind = cost_function_ind + 1;
        ind_pivotToUse = cost_functions_used_correlation_array(cost_function_ind);

        param.pivot = [fullPivotSet(1:ind_pivotToUse-1) '1' fullPivotSet(ind_pivotToUse:end)];
        param.pivotTerm = ind_pivotToUse;
        
        lambda_init = zeros(1, length(J_const_num));
        c_init = [c_init_used(1:ind_pivotToUse-1) c_init_used(ind_pivotToUse+1:end)]'; % use the initial guess passed in from the outside
        
        % convert the matrix-shaped initial value into a vector for lsqlin
%         [z_init, param.feature_count, param.c_array_len] = mat2vec_extra(lambda_init, c_init);
        [z_init, param.matVec_struct] = mat2vec_multi(lambda_init, c_init);

%         [J_cost_use_local, J_const_all_local] = regressorMatrixCondition(z_init, J_cost_use, J_const_all, c_const, param);
        
        [J, A, x, b, J_coeff_out] = ioc_cost_func2(z_init, J_cost_use_reduced, J_const_all, c_const, param);
        
        lsqlin_C = A;
        lsqlin_d = -b; % this is to get Ax+b, since we're doing Ax-(-b)
        lsqlin_A = zeros(size(lsqlin_C, 2), size(lsqlin_C, 2));
        lsqlin_A(1:length(c_init), 1:length(c_init)) = diag(-ones(size(c_init))); % the constrains are Ax<b, so A should be neg for the cost weights
        lsqlin_b = zeros(size(z_init));
        
%         lsqlin_A = zeros(size(linConst_C, 2), size(linConst_C, 2));
%         lsqlin_A = diag(-ones(size(linConst_C, 2))); % the constrains are Ax<b, so A should be neg for the cost weights
        
        lsqlin_Aeq = [];
        lsqlin_beq = [];

        % run lsqlin
        Mesoptions = optimset('Algorithm', 'interior-point', ...
            'TolFun',param_opt.TolFun,'TolX',param_opt.TolX,'TolCon',param_opt.TolCon,...
            'display',displaySetting,'MaxIter',param_opt.maxiter,'Diagnostics','off');
        [z_calc,resnorm_lsqlin_orig,residual_lsqlin,EXITFLAG,output_lsqlin] = lsqlin(lsqlin_C,lsqlin_d,lsqlin_A,lsqlin_b,lsqlin_Aeq,lsqlin_beq,[],[],[],Mesoptions);
        
        % revert the state vector back into array form so we can extract the
        % recovered weights
        [lambda_cell, c_use] = vec2mat_multi(z_calc, param.matVec_struct);
        lambda = lambda_cell{1};
        lambda_out{ind_pivotset} = mat2vec(lambda);
        c_out_temp = zeros(size(param.pivot));
        zc_out_temp = zeros(1, size(param.pivot, 2)-1); zc_out_temp_counter = 0;
        for ind_innerset = 1:length(param.pivot) % copy out the weights
            c_out_temp(ind_innerset) = eval(param.pivot{ind_innerset});
            
            if c_out_temp(ind_innerset) < param.recoverThreshold_c
                c_out_temp(ind_innerset) = 0;
            end
        end
        
        % now normalize the weights
        c_normFactor = sum(c_out_temp);
%         c_out_temp = c_out_temp / c_normFactor;
        
        for ind_innerset = 1:length(param.pivot) % copy out the weights
            if ~strcmpi(param.pivot{ind_innerset}, '1')
                zc_out_temp_counter = zc_out_temp_counter + 1;
                zc_out_temp(zc_out_temp_counter) = c_out_temp(ind_innerset);
            end
        end
        
        c_out_pivot_local = zeros(1, len_cost_functions_full);
        c_out_pivot_local(cost_functions_used_correlation) = c_out_temp; % save the original back in full array
        c_out_pivot_original{ind_pivotset} = c_out_pivot_local;
        
        z_calc_zeroed = z_calc;
        z_calc_zeroed(1:length(zc_out_temp)) = zc_out_temp; % remove the negative values
        z_calc_zeroed = z_calc_zeroed / c_normFactor; % remove the negative values
        c_out_temp = c_out_temp / c_normFactor;
        resnorm_lsqlin = norm(lsqlin_C*z_calc_zeroed - (1/c_normFactor)*lsqlin_d)^2;
        
        c_out_pivot_local = zeros(1, len_cost_functions_full);
        c_out_pivot_local(cost_functions_used_correlation) = c_out_temp; % save the normalized back in full array
        c_out_pivot_normalized{ind_pivotset} = c_out_pivot_local;
        
        % if any columns were removed from coeff, repad them so things will
        % line up properly
        
        % calculate the contribution ratio of the recovered 
        [J_out_contrib{ind_pivotset}, J_out_contrib_ratio{ind_pivotset}, J_out_contrib_array{ind_pivotset}] = ...
            calc_direct_cost(c_out_pivot_local, feature_use_costfctCalc, param);
%         [J, A, x, b, J_coeff_out] = ioc_cost_func2(z_calc, J_cost_use, J_const_all, c_const, param);
        
        % now apply the zeroing and remove entries that are below the
        % threshold
%         c_out_pivot_zeroing = c_out_temp;
%         J_out_pivot_zeroing = J_out_contrib_original{ind_pivotset};
%         c_out_pivot_zeroing(abs(c_out_pivot_zeroing) < param.recoverThreshold_c) = 0;
%         c_out_pivot_zeroing(abs(J_out_pivot_zeroing) < param.recoverThreshold_j) = 0;
        % re-align the cost functions by zero'ing the dropped colinear entries
%         c_out_pivot{ind_pivotset} = c_out_pivot_local; % re-align the cost functions by zero'ing the dropped colinear entries
%         [J_out_contrib{ind_pivotset}, J_out_contrib_ratio{ind_pivotset}, J_out_contrib_array{ind_pivotset}] = ...
%             calc_direct_cost(c_out_pivot{ind_pivotset}, feature_use, param);
%         
%         c_zeroed = [c_out_pivot_zeroing(1:ind_pivotToUse-1) c_out_pivot_zeroing(ind_pivotToUse+1:end)]'; % use the initial guess passed in from the outside
%         [z_zeroed, param.matVec_struct] = mat2vec_multi(lambda, c_zeroed);
%         [J, A, x, b, J_coeff_out] = ioc_cost_func2(z_zeroed, J_cost_use, J_const_all, c_const, param);
        
%         resnorm_constrained = norm(A*x+b)/length(b);
%         resnorm_constrained = resnorm_lsqlin;
        
        % calculate the sign of the dot prod between the pivot and all others
        gradDotSign = [];
        for ind_cgen1 = 1:len_cost_functions_used_correlation
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
        
        costInd = setxor(1:length(J_cost_all), param.pivotTerm);
        
        if 0
            h = figure;
            plot(J);
            hold on; plot(residual_lsqlin, 'r');
            title(['Weights: ' num2str(c_init_used) ' (GT) vs ' num2str(c_out_temp) '(lsqlin)']);
            ylabel(['Pivot (' num2str(ind_pivotToUse) '): ' param.cost_function_names(ind_pivotToUse)]);
            xlabel(['Resnorm: ' num2str(norm(J)) ' (GT) vs ' num2str(norm(residual_lsqlin)) ' (lsqlin)']);
            
            saveas(h, fullfile([saveFolder '_iocFig' num2str(ind_pivotset) '_'  datestr(now, 'yyyymmddHHMMSS') '.fig']));
            close(h);
        end
    end
    
    rref_out = rref_out_full;
    eig_vals = eig_vals_full;
    eig_cutoff = eig_cutoff_full;
%     [~, rref_out, eig_vals, eig_cutoff] = checkCorrelationSVD(J_coeff_out, param.corrThreshold);
    
    % augment J_coeff from matrix removal
%     J_cost_use_reduced = {};
%     J_zero = zeros(size(J_coeff_out(1, :)));
    % zero out the pivots that are correlated
    for ind_pivotset2 = 1:len_cost_functions_full
        if isempty(find(cost_functions_used_correlation == ind_pivotset2, 1))
            J_coeff_out = [J_coeff_out(:, 1:ind_pivotset2-1) J_zero J_coeff_out(:, ind_pivotset2:end)];
        end
    end
    
    % saving the gradient output
    output_inverse(ind_pivotset).J_coeff_individual_out = J_coeff_out;
    output_inverse(ind_pivotset).J_coeff_total_out = J_coeff_total_out;
%     output_inverse(ind_pivotset).resnorm = resnorm_constrained; 
    output_inverse(ind_pivotset).exitflag = EXITFLAG;
    output_inverse(ind_pivotset).A = A;
    output_inverse(ind_pivotset).condA = cond(A);
    output_inverse(ind_pivotset).svd_s = svd(A);
    output_inverse(ind_pivotset).gradDotSign = gradDotSign;
    output_inverse(ind_pivotset).J_out_array = J_out_contrib_array{ind_pivotset};
    output_inverse(ind_pivotset).J_out_contrib = J_out_contrib_ratio{ind_pivotset};
    output_inverse(ind_pivotset).c_original = c_out_pivot_original{ind_pivotset};
    output_inverse(ind_pivotset).c_recovered = c_out_pivot_normalized{ind_pivotset};
    output_inverse(ind_pivotset).c_offset = c_out_pivot_normalized{ind_pivotset}*c_gen_mod(ind_pivotset);
    output_inverse(ind_pivotset).rmse = 0; % this will be updated outside the fct
    output_inverse(ind_pivotset).lambda = lambda_out{ind_pivotset};
    output_inverse(ind_pivotset).c_scalingFactor = c_normFactor;
    output_inverse(ind_pivotset).resnorm_orig = resnorm_lsqlin_orig;
    output_inverse(ind_pivotset).resnorm = resnorm_lsqlin;
    output_inverse(ind_pivotset).resnorm_lsqlin = resnorm_lsqlin;
    output_inverse(ind_pivotset).residual_lsqlin = residual_lsqlin;
    
    output_inverse(ind_pivotset).cost_functions_used_correlation = cost_functions_used_correlation;
    output_inverse(ind_pivotset).rref_out_full = rref_out_full;
    output_inverse(ind_pivotset).eig_vals_full = eig_vals_full;
    output_inverse(ind_pivotset).eig_cutoff_full = eig_cutoff_full;
    output_inverse(ind_pivotset).rref_out = rref_out;
    output_inverse(ind_pivotset).eig_vals = eig_vals;
    output_inverse(ind_pivotset).eig_cutoff = eig_cutoff;
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