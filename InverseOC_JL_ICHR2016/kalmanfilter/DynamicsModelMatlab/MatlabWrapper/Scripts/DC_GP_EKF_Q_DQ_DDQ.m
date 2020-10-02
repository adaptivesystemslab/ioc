classdef DC_GP_EKF_Q_DQ_DDQ < EKF
    %Dual chain constrained gain projection EKF
    
    properties     
        %Constraint Bounds
        upper_state_bound = [];
        lower_state_bound = [];
        upper_mes_bound = [];
        lower_mes_bound = [];
        
       %Model Handle for the second chain
       model_handle = [];
       model_handle2 = [];
       
       D = [];
    end
    
    methods
        function obj = DC_GP_EKF_Q_DQ_DDQ(model,model2)
           obj = obj@EKF;
           obj.model_handle = model;
           obj.model_handle2 = model2;
           
           %Our observations are the sensors attached to the models
            z = [vertcat(model.sensors.measurement); vertcat(model2.sensors.measurement);];
            obj.observation_noise = eye(numel(z));
            %Our state is the q dq ddq of both models
            obj.state = [...
                model.position ; ...
                model.velocity ; ...
                model.acceleration ; ...
                model2.position;
                model2.velocity;
                model2.acceleration ;
                ];
            obj.process_noise = eye(numel(obj.state));
            obj.covariance = eye(numel(obj.state));
            obj.inov_covariance = zeros(numel(z),numel(z));
            obj.sizeZ = numel(z);
            obj.sizeX = numel(obj.state);
            
            %Init Bounds
           obj.upper_state_bound = ones(obj.sizeX,1).*inf;
           obj.lower_state_bound = ones(obj.sizeX,1).*-inf;
           
           obj.upper_mes_bound = ones(obj.sizeZ,1).*inf;
           obj.lower_mes_bound = ones(obj.sizeZ,1).*-inf;
            
        end
        
%         function run_iteration(obj,u,z,timestamp)
%             obj.run_iteration@EKF(u,z);
%             
%             [J1, Jd1] = computeJacobians(obj.model_handle,'frame_ee');
%             [J2, Jd2] = computeJacobians(obj.model_handle2,'frame_ee');
%             obj.D = [zeros(6,size(J1,2)) J1 zeros(6,size(J1,2)) ...
%                 zeros(6,size(J2,2)) -J2 zeros(6,size(J2,2));...
%                 zeros(6,size(J1,2)) Jd1 J1 ...
%                 zeros(6,size(J2,2)) -Jd2 -J2];
%             
%             obj.state=obj.state - obj.D'/(obj.D*obj.D'+1e-12*eye(size(obj.D,1)))*obj.D*obj.state;
%         end
        
        function H = makeH(obj,x)
            obj.model_handle.calculateSensorJacobians();
            H1 = vertcat(obj.model_handle.sensors.obsJacobian);
            obj.model_handle2.calculateSensorJacobians();
            H2 = vertcat(obj.model_handle2.sensors.obsJacobian);
            H = blkdiag(H1,H2);
        end
          
        function x_new = stateUpdate(obj,x,u)
           %State Update
           
            %Update first model
            x_new = x;
            dof = numel(obj.model_handle.joints);
            x_new(1:dof) = x_new(1:dof)+...
                x_new(dof+1:dof*2)*u + x_new(dof*2+1:dof*3)*u^2/2;
            x_new(dof+1:dof*2) = x_new(dof+1:dof*2) + x_new(dof*2+1:dof*3)*u;
            obj.model_handle.position = x_new(1:dof);
            obj.model_handle.velocity = x_new(dof+1:dof*2);
            obj.model_handle.acceleration = x_new(dof*2+1:dof*3);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
            obj.model_handle.forwardAcceleration();
            
            
            
            %Update second model
            dof = numel(obj.model_handle2.joints);            
            b = numel(x_new)-dof*3;
            %Update robot state by integrating
            x_new(b+1:b+dof) = x_new(b+1:b+dof)+...
                x_new(b+dof+1:b+dof*2)*u + x_new(b+dof*2+1:b+dof*3)*u^2/2;
            x_new(b+dof+1:b+dof*2) = x_new(b+dof+1:b+dof*2) + x_new(b+dof*2+1:b+dof*3)*u;
            
            %Forward Kinematics
            obj.model_handle2.position = x_new(b+1:b+dof);
            obj.model_handle2.velocity = x_new(b+dof+1:b+dof*2);
            obj.model_handle2.acceleration = x_new(b+dof*2+1:b+dof*3);
            obj.model_handle2.forwardPosition();
            obj.model_handle2.forwardVelocity();
            obj.model_handle2.forwardAcceleration();
        end
        
        function z = makeMeasure(obj,x)
            %Predict measurement from state
            dof = numel(obj.model_handle.joints);
            obj.model_handle.position = x(1:dof);
            obj.model_handle.velocity = x(dof+1:dof*2);
            obj.model_handle.acceleration = x(dof*2+1:dof*3);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
            obj.model_handle.forwardAcceleration();
            
            dof2 = numel(obj.model_handle2.joints);
            obj.model_handle2.position = x(dof*3+1:dof*3+1+dof2);
            obj.model_handle2.velocity = x(dof*3+dof2+1:dof*3+dof2*2);
            obj.model_handle2.acceleration = x(dof*3+dof2*2+1:dof*3+dof2*3);
            obj.model_handle2.forwardPosition();
            obj.model_handle2.forwardVelocity();
            obj.model_handle2.forwardAcceleration();
            
            %Generate measurements
            z = [vertcat(obj.model_handle.sensors.measurement); vertcat(obj.model_handle2.sensors.measurement)];
        end
        
        function A = makeA(obj,x,u)
            dof1 = numel(obj.model_handle.joints);
            A1 = eye(dof1*3,dof1*3); 
            A1(1:dof1*2,dof1+1:dof1*3) = A1(1:dof1*2,dof1+1:dof1*3)+eye(dof1*2,dof1*2)*u;
            A1(1:dof1,dof1*2+1:dof1*3) = A1(1:dof1,dof1*2+1:dof1*3)+eye(dof1,dof1)*u(1)*u(1)/2;
            
            dof2 = numel(obj.model_handle2.joints);
            A2 = eye(dof2*3,dof2*3); 
            A2(1:dof2*2,dof2+1:dof2*3) = A2(1:dof2*2,dof2+1:dof2*3)+eye(dof2*2,dof2*2)*u;
            A2(1:dof2,dof2*2+1:dof2*3) = A2(1:dof2,dof2*2+1:dof2*3)+eye(dof2,dof2)*u(1)*u(1)/2;
            
            A = blkdiag(A1,A2);
        end
    end
    
    methods(Access=protected)
       
        function K = makeK(obj)
            %Call regular make K
            K = obj.makeK@EKF();
            
            
            %Predict
            x = obj.x_predict+K*obj.dz;
            
            obj.makeMeasure(x);
            %obj.makeH(x);
            
            %Make D matrix
            [J1, Jd1] = computeJacobians(obj.model_handle,'frame_ee');
            [J2, Jd2] = computeJacobians(obj.model_handle2,'frame_ee');
                        
            %J1 = obj.model_handle.sensors(1).baseJacobian;
            %J2 = obj.model_handle2.sensors(1).baseJacobian;
            %obj.D = [zeros(6,size(J1,2)) J1 zeros(6,size(J1,2)) ...
            %    zeros(6,size(J2,2)) -J2 zeros(6,size(J2,2));...
            %    zeros(6,size(J1,2)) Jd1 J1 ...
            %    zeros(6,size(J2,2)) -Jd2 -J2];
            %obj.D = [zeros(6,size(J1,2)) J1 zeros(6,size(J1,2)) ...
            %    zeros(6,size(J2,2)) -J2 zeros(6,size(J2,2))];
            
            %Linearized Constraint on Q
            f1 = obj.model_handle.getFrameByName('frame_ee');
            f2 = obj.model_handle2.getFrameByName('frame_ee');
            T1 = f1.t;
            t1 = T1(1:3,4);
            T2 = f2.t;
            t2 = T2(1:3,4);
            [yd, pd, rd] = dcm2angle((T2(1:3,1:3)*T1(1:3,1:3)')');
            rd = [yd;pd;rd];
            
            %So the constraints are ... 
            obj.D = [J1 zeros(6,size(J1,2)) zeros(6,size(J1,2)) ...
                -J2 zeros(6,size(J2,2)) zeros(6,size(J2,2))];
            b = [t2-t1;rd] + J1*obj.model_handle.position - J2*obj.model_handle2.position;

%             Add Acceleration Constraints
%             obj.D = [obj.D ; [zeros(6,size(J1,2)) Jd1 J1 ...
%                  zeros(6,size(J2,2)) -Jd2 -J2]];
%              b = [b; zeros(6,1)];
%             Add Velocity and Acceleration Constraints             
%              obj.D = [obj.D ; [zeros(6,size(J1,2)) J1 zeros(6,size(J1,2)) ...
%                  zeros(6,size(J2,2)) -J2 zeros(6,size(J2,2));...
%                  zeros(6,size(J1,2)) Jd1 J1 ...
%                  zeros(6,size(J2,2)) -Jd2 -J2]];
%              b = [b; zeros(12,1)];
%             
            
            %Project Gain
            deltaK = obj.D'/(obj.D*obj.D'+1e-6*eye(size(obj.D,1)))*(obj.D*x-b)/(obj.dz'/obj.inov_covariance*obj.dz)*obj.dz'/obj.inov_covariance; 
            %NdK = norm(deltaK)
            %if norm(deltaK) > 0.02
            %    deltaK = deltaK/NdK/50;
            %end
            K = K - deltaK;
            %K = K - D'*pinv(D*D'+1e-8*eye(size(D,1)))*(D*x)/(obj.dz'/obj.inov_covariance*obj.dz)*obj.dz'/obj.inov_covariance;
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