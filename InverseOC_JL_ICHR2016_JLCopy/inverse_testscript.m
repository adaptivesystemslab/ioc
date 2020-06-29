ind_pivotset = 1;

maxiter = 1:30
tol = 1e-16;
lala = [];
lele = [];
for ind_iter = maxiter
param_opt.maxiter = ind_iter;
        param_opt.TolFun = tol;
        param_opt.TolX = tol;
        param_opt.TolCon = tol;
        
        % run lsqlin
        Mesoptions = optimset('Algorithm', 'interior-point', ...
            'TolFun',param_opt.TolFun,'TolX',param_opt.TolX,'TolCon',param_opt.TolCon,...
            'display',displaySetting,'MaxIter',param_opt.maxiter,'Diagnostics','off');
        [z_calc,resnorm_lsqlin_orig,residual_lsqlin,EXITFLAG,output_lsqlin] = lsqlin(lsqlin_C,lsqlin_d,lsqlin_A,lsqlin_b,lsqlin_Aeq,lsqlin_beq,[],[],z_init,Mesoptions);
        
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
        
        c_normFactor = sum(c_out_temp); % now normalize the weights
        
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
        
        lala = [lala; c_out_pivot_local];
        lele = [lele; resnorm_lsqlin];
end