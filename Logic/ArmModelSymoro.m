classdef ArmModelSymoro < GeneralArmModel
    %ARM4DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        % base velocity and acceleration
        baseLinearVelocity = zeros([1,3]); %V0
        baseAngularVelocity = zeros([1,3]); %W0
        baseLinearAcceleration = zeros([1,3]); %VP0
        baseAngularAcceleration = zeros([1,3]); %WP0
        
        %gravity 
        gravity = [0, 0, 0];
    end
    
    methods
        function obj = ArmModelSymoro(jointParams, gravity)       
            obj = obj@GeneralArmModel(jointParams);
            obj.gravity = gravity;
            
        end
                
        function ddq = forwardDynamics(obj, torques, doUpdate)       
            % Function has to be called only after calling update state
            
            if length(obj.joints) == 2
                ddq = forwardDynamics2Dof(obj, torques);
            elseif length(obj.joints) == 3
                ddq = forwardDynamics3Dof(obj, torques);
            else
                ddq = forwardDynamics4Dof(obj, torques);
            end
            
            if ~exist('doUpdate', 'var')
                doUpdate = 0;
            end
               
            if doUpdate
                obj.updateAccelerations(ddq);
            end
            
        end
               
        function tau = inverseDynamics(obj, accelerations)
            
            if length(obj.joints) == 2
                tau = inverseDynamics2Dof(obj, accelerations);
            elseif length(obj.joints) == 3
                tau = inverseDynamics3Dof(obj, accelerations);
            else
                tau = inverseDynamics4Dof(obj, accelerations);
            end
        end
               
        function matrix = getInertialMatrix(obj)
            if length(obj.joints) == 2
                matrix = inertiaMatrix2Dof(obj);
            elseif length(obj.joints) == 3
                matrix = inertiaMatrix3Dof(obj);
            else
                matrix = inertiaMatrix4Dof(obj);
            end
        end
        
        function qVec = getState(obj)
            qVec = zeros(2, obj.totalActuatedJoints);
            for i = 1:obj.totalActuatedJoints
                qVec(1,i) = obj.joints(i).q + obj.joints(i).angleOffset;
                qVec(2,i) = obj.joints(i).dq;
            end
        end
        
        function qVec = getStateWithOffset(obj)
            qVec = zeros(2, obj.totalActuatedJoints);
            for i = 1:obj.totalActuatedJoints
                qVec(1,i) = obj.joints(i).q;
                qVec(2,i) = obj.joints(i).dq;
            end
        end        
        
        function x = getEndEffectorPosition(obj)
            x = [0 0 0];
        end
        
        function result = getDerivativesRelativeToState(obj, control)
            % Computes finite difference derivatives of the forward
            % dynamic equations wrt each state variable
            
            % Dertivatives are computed numerically using central finite
            % differences
                        
            delta = 5e-2;
            % dofs = length(obj.joints);
            dofs = obj.totalActuatedJoints;
            result = zeros(dofs,dofs*2);
            
            currentState = obj.getStateWithOffset();
                        
            for i = 0:dofs-1
                
                changes = zeros(1,dofs);
                changes(i+1) = delta;
                                       
                % Define right and left states for joint i
                leftStateQ = currentState(1,:) - changes;
                rightStateQ = currentState(1,:) + changes;
                leftStateDQ = currentState(2,:) - changes;
                rightStateDQ = currentState(2,:) + changes;

                % Calculate finite difference wrt angle state variable
                obj.updateState(rightStateQ, currentState(2,:));
                ddq = obj.forwardDynamics(control);
                result(:, 2*i+1) = ddq';
                    
                obj.updateState(leftStateQ, currentState(2,:));
                ddq = obj.forwardDynamics(control);
                result(:, 2*i+1) = (result(:, 2*i+1) - ddq')/(2*delta);
                    
                % Calculate finite difference wrt velocity state
                % variable
                obj.updateState(currentState(1,:), rightStateDQ);
                ddq = obj.forwardDynamics(control);
                result(:, 2*i+2) = ddq';
                    
                obj.updateState(currentState(1,:), leftStateDQ);
                ddq = obj.forwardDynamics(control);
                result(:, 2*i+2) = (result(:, 2*i+2) - ddq')/(2*delta);
                  
            end
            obj.updateState(currentState(1,:), currentState(2,:));
        end
    end
end

