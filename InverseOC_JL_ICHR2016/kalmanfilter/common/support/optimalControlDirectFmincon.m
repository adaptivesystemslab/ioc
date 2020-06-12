classdef optimalControlDirectFmincon < optimalControlAbstract % < rlModel
    properties
        docOptStruct;
    end
    
    methods
        function obj = calcOptimalTrajectory(obj, c_cost_doc, c_const_doc)
             obj.modifyModel(); % TODO DEBUG REMOVE LATER           
        
            % set the cost coefficeints on the cost function and constraints
            obj.c_cost = c_cost_doc;
            obj.c_const = c_const_doc;

%             % generate an initial trajectory based on q/dq/ddq constraints
%             obj.featureSet.splineModelInit();
%             obj.featureSet.splineModelSetViaPointsStruct(obj.constPoints);
%             obj.featureSet.calcQFromSpline(obj.windowTime);
            
            % during the setup, featureSet should be prepopulated with the
            % init estimate of the traj
            q_init = obj.featureSet.q;
            
            % update spline model so the time array isn't using the initial
            % trajectory/time, but rather the spline trajectory/time
%             obj.setupSplineViaPoints(obj.featureSet);
%             obj.featureSet.splineModelSetViaPointsStruct(obj.splinePoints);
            
            % extract spline points from the initial trajectory
            y_opt_init_mat = obj.featureSet.splineModel.downloadViaPoints();
            [y_opt_init_vec, obj.matVecStructDoc] = mat2vec_multi(y_opt_init_mat);

            % fmincon setup
            cost_func_handle =  @(y_opt) obj.docCostFunction(y_opt);
            const_func_handle = @(y_opt) obj.docConstFunction(y_opt);
            
            obj.yVecUpdateFeatureSet(y_opt_init_vec);
            [J_cost_sum_preOpt, J_cost_array_preOpt] = calcCostFunction(obj, obj.c_cost, obj.c_const);
            [ceq_vec_preOpt, ceq_q_preOpt, ceq_dq_preOpt, ceq_ddq_preOpt] = calcEqualityConstFunction(obj, obj.c_cost, obj.c_const);

            Mesoptions = optimset('Algorithm','interior-point', 'Jacobian','on', 'FunValCheck','on',...
                'MaxFunEvals',1e8, 'MaxIter',obj.fminconMaxIter, 'Display',obj.fminconDisplaySetting, ... % iter-detailed
                'TolFun',obj.fminconTolFun,'TolX',obj.fminconTolX,'TolCon',obj.fminconTolCon);
            
            % run fmincon
            % actual solver
            [y_out_vec, cost_J, exitFlag, output, lambda, grad, hessian] = fmincon(cost_func_handle, y_opt_init_vec,...
                [],[],[],[],[],[], const_func_handle, Mesoptions);
            
            obj.yVecUpdateFeatureSet(y_out_vec);
            
            [J_cost_sum_postOpt, J_cost_array_postOpt] = calcCostFunction(obj, obj.c_cost, obj.c_const);
            [ceq_vec_postOpt, ceq_q_postOpt, ceq_dq_postOpt, ceq_ddq_postOpt] = calcEqualityConstFunction(obj, obj.c_cost, obj.c_const);

            obj.docOptStruct.y = y_out_vec;
            obj.docOptStruct.cost = cost_J;
            obj.docOptStruct.exitFlag = exitFlag;
            obj.docOptStruct.output = output;
            obj.docOptStruct.lambda = lambda;
            obj.docOptStruct.grad = grad;
            obj.docOptStruct.hessian = hessian;
            
            obj.docOptStruct.c = obj.c_cost;
            obj.docOptStruct.resnorm = norm(ceq_vec_postOpt);
            
            obj.docOptStruct.cost_array = J_cost_array_postOpt;
            obj.docOptStruct.ceq_x_array = ceq_q_postOpt;
            obj.docOptStruct.ceq_dx_array = ceq_dq_postOpt;
            obj.docOptStruct.ceq_ddx_array = ceq_ddq_postOpt;

%             J_cost2 = obj.docConstFunction(y_out_vec);
%             J_const2 = obj.docConstFunction(y_out_vec);

            q_opt = obj.featureSet.q;
            
            if 0 
                h = figure;
                plot(obj.splineTime, q_init, 'r');
                hold on;
                plot(obj.splineTime, q_opt, 'b');
                plot(obj.constPoints.time, obj.constPoints.q, 'x');
                title('Trajectory, quintic spline init (r) vs opt (b)');
                ylabel('Angle [rad]');
                xlabel('Time [s]');
                saveas(h, 'D:\results\ioc\dev\doc.fig');
            end
        end

        % cost function
        function J = docCostFunction(obj, y_opt_vec)
            % update feature set given current opt values
            obj.yVecUpdateFeatureSet(y_opt_vec);
            
            % cost functions
            J = obj.calcCostFunction(obj.c_cost, obj.c_const);
        end
        
        % constrant function
        function [c, ceq] = docConstFunction(obj, y_opt_vec)
            % update feature set given current opt values
            obj.yVecUpdateFeatureSet(y_opt_vec);
            
            % no inequality constraints
            c = [];
            
            % equality constrains are postural
            ceq = obj.calcEqualityConstFunction(obj.c_cost, obj.c_const);
        end
    end
end