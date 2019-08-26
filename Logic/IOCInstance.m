classdef IOCInstance
    %IOCINSTANCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Arm model as defined by Symoro
        dynamicModel;
        % List of names of candidate features
        features;
    end
    
    methods
        function obj = IOCInstance(model, listFeatures)
            %IOCINSTANCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.dynamicModel = model;
            obj.features = listFeatures;
        end
        
        function df_dx = getDynamicsStateDerivatives(obj, dt, control)
            % This function computes Fx = df/dx
            % Derivates are computed with respect to state-space
            % representation, which is defined as follows:
            %   x1 = q, x2 = dq
            %   dx1 = x2 = dq, dx2 = ddq, with ddq given by M(q)^-1(torque - H(q,dq))
            %   
            % |dx1| = |x2|
            % |dx2|   |ddq|
            
            % In discrete representation
            % |dx1| = |q_t + dt*dq_t  |
            % |dx2|   |dq_t + dt*qdq_t|
            
            % Define final matrix. For each dof, there are two state
            % equations
            dofs = obj.dynamicModel.totalActuatedJoints;
            df_dx = eye(dofs*2, dofs*2);
            
            ddqDerivatives = obj.dynamicModel.getDerivativesRelativeToState(control);
            ddqDerivatives = dt * ddqDerivatives;
            
            for i = 0:dofs-1
                row_q = 2*i+1;
                row_dq = 2*i+2;
                df_dx(row_q, row_dq) = dt;
                df_dx(row_dq,:) = df_dx(row_dq,:) + ddqDerivatives(i+1,:);
            end
            
        end
        
        function df_du = getDynamicsControlDerivatives(obj, dt, control)
            
            inertiaMatrix = obj.dynamicModel.getInertialMatrix();
            inv = pinv(inertiaMatrix);
            
            % Define final matrix. For each dof, there are two state
            % equations
            dofs = obj.dynamicModel.totalActuatedJoints;
            df_du = zeros(dofs*2, dofs);
            
            for i = 0:dofs-1
                row_q = 2*i+1;
                row_dq = 2*i+2;
                
                df_du(row_q,:) = zeros(1,dofs);
                u = zeros(1,dofs)';
                u(i+1) = 1;
                df_du(row_dq,:) = (dt*inv*u)';
                % df_du(2:2:end,i+1) = (dt*inv*u)';
                
            end
                        
        end
                    
        function dp_dx = getFeaturesStateDerivatives(obj, velocities)
            
            % In the refactorized implementation, this function will
            % consider a list of all features to be computed
            
            % For the moment, these derivatives are hard coded and we are
            % only considering torque and velocities as features
            
            % CODE BELOW: To be used when velocities are included as cost
            % function
%             dp_dq = zeros(8,4); % derivatives wrt to angles q
%             dp_dv = eye(4,4) * (2*velocities)'; % derivatives wrt velocities dq
%             dp_dv = [zeros(1,4); dp_dv(1,:); zeros(1,4); dp_dv(2,:);...
%                 zeros(1,4); dp_dv(3,:); zeros(1,4); dp_dv(4,:)];
%             
%             dp_dx_next = [dp_dq dp_dv];
            
            % Current trajectory only has torque as cost function, thus
            % these derivatives are all zero
            dofs = obj.dynamicModel.totalActuatedJoints;
            dp_dx = zeros(dofs,dofs*2);
            
        end
        
        function matrix = getFeaturesControlDerivatives(obj, control)
            dofs = obj.dynamicModel.totalActuatedJoints;
  
            dp1_du = zeros(dofs,dofs); %Velocities features' derivatives wrt control
            dp2_du = eye(dofs,dofs);
            dp2_du(1:dofs+1:end) = 2*control;
            
            % Use when velocities are considered as cost function.
            % matrix = [dp1_du dp2_du];            
            matrix = dp2_du;            
        end
        
        
        function [df_dx, df_du, dp_dx, dp_du] = getDerivativesNewObservation(obj,...
            control, angles, velocities, dt)
            
            % First, update model with new trajectory point
            obj.dynamicModel.updateState(angles, velocities);
            % Get derivatives of dynamics wrt to state vector
            df_dx = obj.getDynamicsStateDerivatives(dt, control);
            
            % First, update model with new trajectory point
            obj.dynamicModel.updateState(angles, velocities);
            % Get derivatives of dynamics wrt to control signal
            df_du = obj.getDynamicsControlDerivatives(dt, control);
            
            % First, update model with new trajectory point
            obj.dynamicModel.updateState(angles, velocities);
            % Get derivatives of cost functions wrt to state vector
            dp_dx = obj.getFeaturesStateDerivatives(velocities);
            
            % First, update model with new trajectory point
            obj.dynamicModel.updateState(angles, velocities);
            % Get derivatives of cost functions wrt to control signal
            dp_du = obj.getFeaturesControlDerivatives(control);            
        end
    end
end

