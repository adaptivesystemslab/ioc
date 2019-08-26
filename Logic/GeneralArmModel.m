classdef GeneralArmModel < handle
    % This class contains all kinematics and dynamics parameters required to define 
    % an articulated chain. It also includes forward and inverse
    % dynamic and kinematic functions
    
    properties
        joints;
        totalActuatedJoints=0;
    end
    
    methods(Abstract=true)
        % Functions be implemented by the subclassess
        ddq = forwardDynamics(obj, torques, doUpdate);
        tau = inverseDynamics(obj, accelerations);
        matrix = getInertialMatrix(obj);
        x = getEndEffectorPosition(obj);
        result = getDerivativesRelativeToState(obj, control)
    end     
    
    methods
        function obj = GeneralArmModel(jointParams)
            obj.joints = ModelJoint(jointParams(1));
            for i = 2:length(jointParams)
                obj.joints(i) = ModelJoint(jointParams(i));
            end
            
            for i = 1:length(jointParams)
                if ~obj.joints(i).endEffector
                    obj.totalActuatedJoints = obj.totalActuatedJoints+1;
                end
            end
        end
        
        function updateState(obj, newQ, newDQ)
            if ~(isempty(newQ) && isempty(newDQ))
                for i= 1:length(newQ)
                    obj.joints(i).updateState(newQ(i), newDQ(i));
                end
            end
        end
                
        function qVec = getState(obj)
            qVec = zeros(2, length(obj.joints));
            for i = 1:length(obj.joints)
                qVec(1,i) = obj.joints(i).q;
                qVec(2,i) = obj.joints(i).dq;
            end
        end
        
        function updateAccelerations(obj, ddq)            
            for i = 1:obj.totalActuatedJoints
                obj.joints(i).ddq = ddq(i);
            end
        end
        
        function ddqVec = getAccelerations(obj)        
            ddqVec = zeros(1, obj.totalActuatedJoints);
            for i = 1:obj.totalActuatedJoints
                ddqVec(i) = obj.joints(i).ddq;
            end
        end                
    end
   
end

