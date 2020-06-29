classdef MC_EKF_Q_DQ_DDQ < EKF
    % Multichain EKF with state projection constraints using
    
    properties
        %List of model handles, we can add models to the EKF which will
        %update the size of state, cov, and mes, and ect
        model_handles = rlCModel.empty();
        %Keeps track of where the state starts for each model
        model_state_indices = 1;
        
        %Constraints on the models
        constraints = EKF_Constraint.empty();
        
        %Initial covariance for q dq and ddq
        icov = [1 1 1];
        
        %Unprojected state
        state_orig;
        
        
        %Process noise for accleration
        eta_r = 500000;
        eta_p = 500000;
        dt = 1/100;
    end
    
    methods
        function addModel(obj,model)
            %Add a model to the EKF
            
            %Add model to the model list
            obj.model_handles(end+1) = model;
            
            %Pull out prismatic vs revolute joints
            p_joint_indxs = contains({model.joints.type},'prismatic');
            r_joint_indxs = contains({model.joints.type},'revolute');
            
            %Update state
            obj.state = [obj.state; model.position; model.velocity; model.acceleration];
            
            %Update covariance
            obj.covariance = blkdiag(obj.covariance,eye(numel(model.position)*3));
            
            %Update process noise
            %Process noise of EKF left at identity
            %This eta is result of the optimization
            dim = numel(model.position);
            G = [];
            for i=3:-1:1
                G = [G; ones(dim,1)*obj.dt^i / factorial(i)];
            end
            P_tmp = G*G';
            P = zeros(size(P_tmp));
            for i=1:3
                for j=i:3
                    D = diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim));
                    D(p_joint_indxs) = D(p_joint_indxs)*obj.eta_p;
                    D(r_joint_indxs) = D(r_joint_indxs)*obj.eta_r;
                    P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
                        diag(D);
                end
            end
            P = P+P' - diag(diag(P));
            %P(1:dim,1:dim) = P(1:dim,1:dim) +eye(dim)*0.00001;
            %P(dim+1:dim*2,dim+1:dim*2) = P(dim+1:dim*2,dim+1:dim*2) +eye(dim)*0.01;
            %P(dim*2+1:end,dim*2+1:end) = P(dim*2+1:end,dim*2+1:end) +eye(dim)*0.1;
            obj.process_noise = blkdiag(obj.process_noise,P);
            
            z = vertcat(model.sensors.measurement);
            obj.sizeZ = obj.sizeZ + numel(z);
            obj.sizeX = obj.sizeX + numel(model.position)*3;
            
            %Add the starting point of next model state to index
            obj.model_state_indices(end+1) = obj.model_state_indices(end)+...
                numel(model.position)*3;
            
            %Default observation noise
            obj.observation_noise = blkdiag(obj.observation_noise,eye(numel(z)));
        end
        
        function run_iteration(obj,u,z,timestamp)
            %Run ekf iteration and project state onto constrained space
            
            %Run the regular EKF update
            run_iteration@EKF(obj,u,z);
            
            %Build the constraints and project state gupta 2007
            
            %Build A and b for position constraints
            num_constraints = numel(obj.constraints);
            A_cell = cell(num_constraints,1);
            b_cell = cell(num_constraints,1);
            for i=1:num_constraints
                %If it is a Jacobian based constraint we have to figure out
                %what 2 models to set it for
                C = obj.constraints(i);
                %Find the starting index of each model
                dof1 = numel(C.models(1).position);
                indx1 = obj.model_state_indices(find(C.models(1) == obj.model_handles,1));
                dof2 = numel(C.models(2).position);
                indx2 = obj.model_state_indices(find(C.models(2) == obj.model_handles,1));
                if(isa(C,'EKF_Constraint_BF') || isa(C,'EKF_Constraint_T0EE') || isa(C,'EKF_Constraint_BF_vel'))
                    
                    %Shift index for velocity constraints
                    if isa(C,'EKF_Constraint_BF_vel')
                       indx1 = indx1+dof1;
                       indx2 = indx2+dof2;
                    end
                    
                    %Build the A for this particular constraint
                    A_c = C.A;
                    A = zeros(size(A_c,1),obj.sizeX);
                    if indx1 == indx2
                        %The constraint is within the same model
                        A = zeros(size(A_c,1),obj.sizeX);
                        A(:,indx1:indx1+dof1-1) = A_c;
                    else
                        %The constraint is between 2 different models
                        A(:,indx1:indx1+dof1-1) = A_c(:,1:dof1);
                        A(:,indx2:indx2+dof2-1) = A_c(:,dof1+1:end);
                    end
                    
                    %For stacking
                    A_cell{i} = A;
                    b_cell{i} = C.b;
                else
                    
                    
                end
            end
            A_pos = vertcat(A_cell{:});
            b_pos = vertcat(b_cell{:});
            
            if(~isempty(A_pos))
                W_k_inv  = obj.covariance;
                Y = W_k_inv*A_pos'/(A_pos*W_k_inv*A_pos' + eye(size(A_pos,1))*5e-7);
                %Save unprojected state
                obj.state_orig = obj.state;
                %Project State
                obj.state = obj.state - Y*(A_pos*obj.state - b_pos);
                %Update covariance
                obj.covariance = (eye(size(obj.covariance)) - Y*A_pos)*obj.covariance;
                %Maybe we don't need this
                obj.covariance = nearestSPD(obj.covariance);
            end
        end
        
        function H = makeH(obj,x)
            %Build Jacobian
            
            %Calcualte all sensor Jacobians for each model and concatenate
            %the sensor Jacobians
            num_mdls = numel(obj.model_handles);
            H_cell = cell(num_mdls,1);
            for i=1:num_mdls
                obj.model_handles(i).calculateSensorJacobians();
                H_cell{i} = vertcat(obj.model_handles(i).sensors.obsJacobian);
            end
            %Concat diag
            H = blkdiag(H_cell{:});
        end
        
        function x_new = stateUpdate(obj,x,u)
            %State Update
            x_new = x;
            %Update models
            for i=1:numel(obj.model_handles)
                dof = numel(obj.model_handles(i).joints);
                indx = obj.model_state_indices(i);
                x_new(indx:indx+dof-1) = x_new(indx:indx+dof-1)+...
                    x_new(indx+dof:indx+dof*2-1)*u + x_new(indx+dof*2:indx+dof*3-1)*u^2/2;
                x_new(indx+dof:indx+dof*2-1) = x_new(indx+dof:indx+dof*2-1) +...
                    x_new(indx+dof*2:indx+dof*3-1)*u;
                obj.model_handles(i).position = x_new(indx:indx+dof-1);
                obj.model_handles(i).velocity = x_new(indx+dof:indx+dof*2-1);
                obj.model_handles(i).acceleration = x_new(indx+dof*2:indx+dof*3-1);
                obj.model_handles(i).forwardPosition();
                obj.model_handles(i).forwardVelocity();
                obj.model_handles(i).forwardAcceleration();
            end
        end
        
        function z = makeMeasure(obj,x)
            %Predict measurement from state
            
            num_models = numel(obj.model_handles);
            z_cell = cell(num_models,1);
            for i=1:num_models
                dof = numel(obj.model_handles(i).joints);
                indx = obj.model_state_indices(i);
                obj.model_handles(i).position = x(indx:indx+dof-1);
                obj.model_handles(i).velocity = x(indx+dof:indx+dof*2-1);
                obj.model_handles(i).acceleration = x(indx+dof*2:indx+dof*3-1);
                obj.model_handles(i).forwardPosition();
                obj.model_handles(i).forwardVelocity();
                obj.model_handles(i).forwardAcceleration();
                z_cell{i} = vertcat(obj.model_handles(i).sensors.measurement);
            end
            z = vertcat(z_cell{:});
        end
        
        function addConstraint(obj,c)
            %Add a constraint to the EKF
            obj.constraints(end+1) = c;
        end
        
        function A = makeA(obj,x,u)
            num_models = numel(obj.model_handles);
            A_cell = cell(num_models,1);
            for i=1:num_models
                dof = numel(obj.model_handles(i).joints);
                A = eye(dof*3,dof*3);
                A(1:dof*2,dof+1:dof*3) = A(1:dof*2,dof+1:dof*3)+eye(dof*2,dof*2)*u;
                A(1:dof,dof*2+1:dof*3) = A(1:dof,dof*2+1:dof*3)+eye(dof,dof)*u(1)*u(1)/2;
                A_cell{i}=A;
            end
            A = blkdiag(A_cell{:});
        end
        
    end
    
    
end