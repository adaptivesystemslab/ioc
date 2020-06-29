classdef optimalControlInverseKKT < optimalControlAbstract % < rlModel
    properties
        h = [0 -1/2 1/2]*1e-6;
        recoverThreshold_c = 1e-12;

        iocOptStruct;
        iocOptNormStruct;
        optPivotInd;
    end
    
    methods
        function obj = calcOptimalWeights(obj)
            % initialize
            c_init = zeros(size(obj.costFunctionStruct));
            
            % calculate gradient for J_cost and J_const
            [J_cost, grad_cost, grad_const] = obj.calcGradient();
            
            grad_cost_check = [grad_cost{:}];
            grad_const_check = [grad_const{:}];
            
            % iterate through the pivots to find the recovered weights for all
            [obj.iocOptStruct, obj.iocOptNormStruct] = obj.pivotGen(c_init, grad_cost, grad_const);
            
            % select the best pivot
            obj.optPivotInd = obj.pivotResolve(obj.iocOptNormStruct);
        end
        
        function [J_cost, grad_cost, grad_const] = calcGradient(obj)
%             obj.featureSet.calcQFromSpline(obj.windowTime);
%             obj.featureSet.calcFeaturesFromQ();
            
            y_opt_mat = obj.featureSet.splineModel.downloadViaPoints();
            [y_opt_vec, obj.matVecStructDoc] = mat2vec_multi(y_opt_mat);
            
            len_ycombo = obj.matVecStructDoc.feature_count*obj.matVecStructDoc.feature_width*obj.matVecStructDoc.entry_count;
            
            J_template = zeros(len_ycombo, 1);
            [~, J_cost_array] = obj.calcCostFunction(ones(size(obj.costFunctionStruct)), obj.c_const);
            [ceq] = obj.calcEqualityConstFunction(obj.c_cost, obj.c_const);
            len_h = numel(obj.h);
            len_cost = numel(J_cost_array);
            len_const = numel(ceq);
            
            J_h_array = cell(1, len_cost);
            J_h_array(:) = {J_template}; % preallocate the cost function array
            g_h_array = cell(1, len_const);
            g_h_array(:) = {J_template};
            grad_cost = cell(1, len_cost);
            grad_cost(:) = {J_template};
            grad_const = cell(1, len_const);
            grad_const(:) = {J_template};
            
%             figure;
%             hold on
            
            for ind_knots = 1:len_ycombo
                for ind_h = 1:length(obj.h)
                    y_opt_vec_offset = y_opt_vec;
                    y_opt_vec_offset(ind_knots) = y_opt_vec_offset(ind_knots) + obj.h(ind_h);
                    
                    % calculate features
                    obj.yVecUpdateFeatureSet(y_opt_vec_offset);
                    
%                     hold on
%                     subplot(311);
%                     plot(obj.featureSet.q);
%                      hold on
%                     subplot(312);
%                     plot(obj.featureSet.dq);
%                      hold on
%                     subplot(313);
%                     plot(obj.featureSet.ddq);
% %                     pause(0.01);
                    
                    % extract the cost/const functions
                    [~, J_cost_array] = obj.calcCostFunction(ones(size(obj.costFunctionStruct)), obj.c_const);
                    [ceq] = obj.calcEqualityConstFunction(obj.c_cost, obj.c_const);
                    
                    for ii_J = 1:len_cost
                        J_h_array{ii_J}(ind_knots, ind_h) = J_cost_array(ii_J);
                    end
                    
                    % and now build the constraint functions
                    for ii_J = 1:len_const
                        g_h_array{ii_J}(ind_knots, ind_h) = ceq(ii_J);
                    end
                end
            end
                
            for ind_knots = 1:len_ycombo
                % iterate through all th entries in the J array, as calculated above,
                % to note how they change corresponding to the changes in q_knots
                for ii_J = 1:len_cost
                    J_curr_knot = J_h_array{ii_J}(ind_knots, :); % this value is normalized since the normalization function is in the cost calculations
                    grad_cost{ii_J}(ind_knots) = obj.calcDiffNum(J_curr_knot);
                end
                
                % need to break the full length array into 3 smaller components to
                % check which (ie x, dx or ddx) is currently active to determine where
                % the constraints are
                for ii_J = 1:len_const
                    g_curr_knot = g_h_array{ii_J}(ind_knots, :); % this value is normalized since the normalization function is in the cost calculations
                    grad_const{ii_J}(ind_knots) = obj.calcDiffNum(g_curr_knot);
                end
            end
            
            % pull out the actual cost function. at h=0, it shouldn't have shifted
            J_cost = zeros(1, len_cost);
            h_ind1 = obj.h == 0;
            for ind = 1:len_cost
                J_cost(ind) = J_h_array{ind}(1, h_ind1);
            end
        end
        
        function J = calcDiffNum(obj, x)
            h_ind1 = obj.h == obj.h(2);
            h_ind2 = obj.h == obj.h(3);
            
            num = x(:, h_ind1) - x(:, h_ind2);
            den = obj.h(h_ind1) - obj.h(h_ind2);
            num_diff = num / den; % sum of (ddq(x+h) - ddq(x-h))/2h
            J = sum(num_diff, 2) / size(x, 1); % summing over all units of time, divided by time
        end
        
        function [iocOptStruct, iocOptNormStruct] = pivotGen(obj, c_init_in, grad_cost, grad_const)
            len_costFunctionsUsed = length(obj.costFunctionStruct);
            
            for ind_pivot = 1:len_costFunctionsUsed
                % init guess before lsqlin
                lambda_init = zeros(length(grad_const), 1);
                c_init = [c_init_in(1:ind_pivot-1) c_init_in(ind_pivot+1:end)]'; 
                
                % convert the matrix-shaped initial value into a vector for lsqlin
                [z_init, obj.matVecStructIoc] = mat2vec_multi(lambda_init, c_init);
                
                % calculate the required matrices for IOC
                [J_init, A_init, x_init, b_init, ~] = obj.ioc_cost_func(z_init, grad_cost, grad_const, obj.c_const, ind_pivot);
                C_lsqlin = A_init;
                d_lsqlin = -b_init; % this is to get Ax+b, since we're doing Ax-(-b)
                A_lsqlin = zeros(size(C_lsqlin, 2), size(C_lsqlin, 2));
                A_lsqlin(1:length(c_init), 1:length(c_init)) = diag(-ones(size(c_init))); % the constrains are Ax<b, so A should be neg for the cost weights
                b_lsqlin = zeros(size(z_init));

                Aeq_lsqlin = [];
                beq_lsqlin = [];
                
                % run IOC via lsqlin
                Mesoptions = optimset('Algorithm','active-set', 'Diagnostics','off',...
                    'MaxFunEvals',1e8, 'MaxIter',obj.fminconMaxIter, 'Display',obj.fminconDisplaySetting, ... % iter-detailed
                    'TolFun',obj.fminconTolFun,'TolX',obj.fminconTolX,'TolCon',obj.fminconTolCon);
  
                [z_lsqlin,resnorm_lsqlin,residual_lsqlin,exitFlag,output] = lsqlin(C_lsqlin,d_lsqlin,A_lsqlin,b_lsqlin,Aeq_lsqlin,beq_lsqlin,[],[],z_init,Mesoptions);
                
                % save the lsqlin output
                c_normFactor = 1;
                iocOptStruct(ind_pivot) = obj.calcLsqLinOutput(c_normFactor, z_lsqlin, C_lsqlin, d_lsqlin, ind_pivot);
                
                % now normalize the calculated values   
                c_normFactor = sum(iocOptStruct(ind_pivot).c);
                iocOptNormStruct(ind_pivot) = obj.calcLsqLinOutput(c_normFactor, z_lsqlin, C_lsqlin, d_lsqlin, ind_pivot);
            end
        end
        
        function [J, A, x, b, J_cost_combined] = ioc_cost_func(obj, z_init, grad_cost, grad_const, c_const, pivotTerm)
            % IOC cost function calc. combine the matrices generated by setup_ioc into
            % the form that the pivot requires
            [lambda_use_cell, c_use] = vec2mat_multi(z_init, obj.matVecStructIoc);
            lambda_use = mat2vec(lambda_use_cell{1});
            
            J_cost_combined = horzcat(grad_cost{:});
            J_const_combined = horzcat(grad_const{:});
            
            % calculate the cost function side of the eq'n
            J_b = grad_cost{pivotTerm}; % pull out the pivot term
            grad_cost{pivotTerm} = []; % the nth term is the pivot, don't include it (ie J_coeff_ddq)
            A_cost = horzcat(grad_cost{:});
            
            % calculate the const function side fo the eq'n
            % J_coeff_constBlk = blkdiag(J_coeff_const{:});
            % A = blkdiag(A_cost, zeros(size(J_coeff_constBlk))); % include the constraints base
            % b = [J_b; zeros(size(J_coeff_constBlk, 1), 1)];
            
            A = A_cost; % have not included the constraints
            b = J_b;
            
            % now deal with the constraints. merge the two matrices
            A_const = horzcat(grad_const{:});
            
             % TODO mess with J_const_combined
            J_const_abs = abs(A_const);
            A_const(J_const_abs < 1e-5) = 0;
            A_const(J_const_abs > 1e-3) = 1;
            
            A = [A A_const];
            
            %  set up the x array and its weights (1 for c, c_const for lambda)
            len_c = 1:length(c_use);
            len_constx = len_c(end) + (1:numel(obj.constPoints.q));
            len_constdx = len_constx(end) + (1:numel(obj.constPoints.dq));
            len_constddx = len_constdx(end) + (1:numel(obj.constPoints.ddq));
            x_weight = zeros(size(z_init));
            x_weight(len_c) = 1;
            x_weight(len_constx) = c_const(1);
            x_weight(len_constdx) = c_const(2);
            x_weight(len_constddx) = c_const(3);
            
            x = x_weight .* [c_use; lambda_use];
            
            J = A*x + b;
            
            % check for length and width of A
            heightX = length(x);
            heightA = size(J_cost_combined, 1);
            
            if heightA < heightX % the A matrix is 'wide'. apply regulariation
                J = J + (1e-9)*norm(x);
            else % else square, or 'tall' matrix
                % don't do anything
            end
        end
        
        function optStruct = calcLsqLinOutput(obj, c_normFactor, z, C, d, ind_pivot)
            % revert the state vector back into array form so we can extract the
            % recovered weights
            [lambda_cell, c_use] = vec2mat_multi(z, obj.matVecStructIoc);
            lambda = mat2vec(lambda_cell{1});
            c = [c_use(1:ind_pivot-1); 1; c_use(ind_pivot:end)]'; % insert the pivot value ('1') back in
            c(c < obj.recoverThreshold_c) = 0; % zero out values that are too small
            
            % apply normalization
            z_norm = z/c_normFactor;
            c_norm = c/c_normFactor;
            d_norm = d/c_normFactor;
            lambda_norm = lambda/c_normFactor;
            
            % calculate the cost array of the recovered
            [J_cost, J_cost_array] = obj.calcCostFunction(c_norm, obj.c_const);
            
            % save to a struct
            optStruct.z = z_norm;
            optStruct.c = c_norm;
            optStruct.lambda = lambda_norm;
            optStruct.resnorm = norm(C*z_norm - d_norm)^2;
            optStruct.J_cost = J_cost;
            optStruct.J_cost_array = J_cost_array;
        end
        
        function optPivotInd = pivotResolve(obj, optNormStruct)
            errorScore = 1e6*ones(size(optNormStruct)); % init to a large number
            t_input = obj.featureSet.time;
            q_input = obj.featureSet.q;
            
            switch obj.pivotMethod
                case 'rmse'
                    % call the direct problem
                    docSim = optimalControlSim();
                    docSim.initOptimalControl(obj.featureSet); % load a init trajectory (splined, in this case)
                    
                    for ind_pivot = 1:length(optNormStruct)
                        fprintf('Solving the direct problem using pivot weights recovered (act: %u/%u)..\n', ind_pivot, length(optNormStruct));
                
                        % regenerating the traj with the recovered weights
                        c_cost = optNormStruct(ind_pivot).c;
                        c_const = obj.c_const;
                        docSim.calcOptimalTrajectory(c_cost, c_const);
                        
                        t_sim = docSim.featureSet.time;
                        q_sim = docSim.featureSet.q;
                        
                        errorScore(ind_pivot) = obj.rmseFct(q_input, q_sim);
                        
                        if 1
                            figure; 
                            plot(t_input, q_input);
                            hold on;
                            plot(t_sim, q_sim, '--');
                            title(['q input (solid) vs q rms (dashed), RMSE = ' num2str(errorScore(ind_pivot))]);
                        end
                    end

                case 'resnorm'
                    errorScore = [optNormStruct.resnorm];
            end
            
            [~, optPivotInd] = min(errorScore);
        end
    end
end