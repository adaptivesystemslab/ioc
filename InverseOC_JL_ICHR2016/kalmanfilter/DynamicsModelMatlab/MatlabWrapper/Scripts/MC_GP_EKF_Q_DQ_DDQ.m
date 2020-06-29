classdef MC_GP_EKF_Q_DQ_DDQ < EKF
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
        
        %Process noise for accleration
        eta = 10;
    end
    
    methods
        function addModel(obj,model)
            %Add a model to the EKF
            
            %Add model to the model list
            obj.model_handles(end+1) = model;
            
            %Update state
            obj.state = [obj.state; model.position; model.velocity; model.acceleration];
            
            %Update covariance
            obj.covariance = blkdiag(obj.covariance,eye(numel(model.position)*3));
            
            %Update process noise
            %Process noise of EKF left at identity
            %This eta is result of the optimization
            dt = 0.01;
            dim = numel(model.joints);
            G = [];
            for i=3:-1:1
                G = [G; ones(dim,1)*dt^i / factorial(i)];
            end
            P_tmp = G*G'*10;
            P = zeros(size(P_tmp));
            for i=1:3
                for j=i:3
                    P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
                        diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
                end
            end
            P = P+P' - diag(diag(P));
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
%             
%             %Build A and b
%             num_constraints = numel(obj.constraints);
%             A_cell = cell(num_constraints,1);
%             b_cell = cell(num_constraints,1);
%             for i=1:num_constraints
%                 %If it is a Jacobian based constraint we have to figure out
%                 %what 2 models to set it for
%                 C = obj.constraints(i);
%                 if(isa(C,'EKF_Constraint_BF') || isa(C,'EKF_Constraint_BFConst'))
%                     
%                     %Find the starting index of each model
%                     dof1 = numel(C.models(1).position);
%                     indx1 = obj.model_state_indices(find(C.models(1) == obj.model_handles,1));
%                     dof2 = numel(C.models(2).position);
%                     indx2 = obj.model_state_indices(find(C.models(2) == obj.model_handles,1));
%                     
%                     %Build the A for this particular constraint
%                     A_c = C.A;
%                     A = zeros(size(A_c,1),obj.sizeX);
%                     A(:,indx1:indx1+dof1-1) = A_c(:,1:dof1);
%                     A(:,indx2:indx2+dof2-1) = A_c(:,dof1+1:end);
%                     A_cell{i} = A;
%                     b_cell{i} = C.b;
%                 else
%                     
%                     
%                 end
%             end
%             A = vertcat(A_cell{:});
%             b = vertcat(b_cell{:});
%             
%             if(~isempty(A))
%                 %For now lets assume weighting matrix is I
%                 %TODO !!!!!! investigate why it should be P_k|k^-1
%                 W_k = eye(size(obj.A));
%                 W_k_inv = eye(size(obj.A));
%                 Y = W_k_inv*A'/(A*W_k_inv*A' + eye(size(A,1))*1e-6);
%                 %Project State
%                 obj.state = obj.state - Y*(A*obj.state - b);
%                 %Update covariance
%                 obj.covariance = (eye(size(obj.covariance)) - Y*A)*obj.covariance;
%                 obj.covariance = nearestSPD(obj.covariance);
%                 
%             end
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
    
    methods(Access=protected)
        
        function K = makeK(obj)
            %Call regular make K
            K = obj.makeK@EKF();
            %Predict
            x = obj.x_predict+K*obj.dz;
            obj.makeMeasure(x);
            
            %Build A and b
            num_constraints = numel(obj.constraints);
            A_cell = cell(num_constraints,1);
            b_cell = cell(num_constraints,1);
            for i=1:num_constraints
                %If it is a Jacobian based constraint we have to figure out
                %what 2 models to set it for
                C = obj.constraints(i);
                if(isa(C,'EKF_Constraint_BF') || isa(C,'EKF_Constraint_BFConst'))
                    
                    %Find the starting index of each model
                    dof1 = numel(C.models(1).position);
                    indx1 = obj.model_state_indices(find(C.models(1) == obj.model_handles,1));
                    dof2 = numel(C.models(2).position);
                    indx2 = obj.model_state_indices(find(C.models(2) == obj.model_handles,1));
                    
                    %Build the A for this particular constraint
                    A_c = C.A;
                    A = zeros(size(A_c,1),obj.sizeX);
                    A(:,indx1:indx1+dof1-1) = A_c(:,1:dof1);
                    A(:,indx2:indx2+dof2-1) = A_c(:,dof1+1:end);
                    A_cell{i} = A;
                    b_cell{i} = C.b;
                else
                    
                    
                end
            end
            A = vertcat(A_cell{:});
            b = vertcat(b_cell{:});
            
            %Project Gain
            deltaK = A'/(A*A'+1e-6*eye(size(A,1)))*(A*x-b)/(obj.dz'/obj.inov_covariance*obj.dz)*obj.dz'/obj.inov_covariance;
            K = K - deltaK;
        end
        
        function P = makeP(obj)
            P = obj.P_predict-obj.K*obj.H*obj.P_predict;
            P = obj.nearestSPD(P);
        end
        
        function Ahat = nearestSPD(obj,A)
            % nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
            % usage: Ahat = nearestSPD(A)
            %
            % From Higham: "The nearest symmetric positive semidefinite matrix in the
            % Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
            % where H is the symmetric polar factor of B=(A + A')/2."
            %
            % http://www.sciencedirect.com/science/article/pii/0024379588902236
            %
            % arguments: (input)
            %  A - square matrix, which will be converted to the nearest Symmetric
            %    Positive Definite Matrix.
            %
            % Arguments: (output)
            %  Ahat - The matrix chosen as the nearest SPD matrix to A.
            
            % test for a square matrix A
            [r,c] = size(A);
            if r ~= c
                error('A must be a square matrix.')
            elseif (r == 1) && (A <= 0)
                % A was scalar and non-positive, so just return eps
                Ahat = eps;
                return
            end
            
            % symmetrize A into B
            B = (A + A')/2;
            
            % Compute the symmetric polar factor of B. Call it H.
            % Clearly H is itself SPD.
            [U,Sigma,V] = svd(B);
            H = V*Sigma*V';
            
            % get Ahat in the above formula
            Ahat = (B+H)/2;
            
            % ensure symmetry
            Ahat = (Ahat + Ahat')/2;
            
            % test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
            p = 1;
            k = 0;
            while p ~= 0
                [R,p] = chol(Ahat);
                k = k + 1;
                if p ~= 0
                    % Ahat failed the chol test. It must have been just a hair off,
                    % due to floating point trash, so it is simplest now just to
                    % tweak by adding a tiny multiple of an identity matrix.
                    mineig = min(eig(Ahat));
                    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
                end
            end
        end
        
    end
    
    
end