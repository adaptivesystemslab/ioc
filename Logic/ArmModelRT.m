classdef ArmModelRT < GeneralArmModel

    properties
        model;
    end
    
    methods(Access=private)
        % create model according to Peter Corke Toolbox
        function instantiateModel(obj)
            L = repmat(RevoluteMDH(),1, length(obj.joints)); 
            for i = 1:length(obj.joints)
                currJoint = obj.joints(i);
                L(i) = RevoluteMDH('d', currJoint.offset, 'a', currJoint.length,...
                    'alpha', currJoint.twist, 'offset', currJoint.angleOffset,...
                    'm', currJoint.mass, 'r', currJoint.com, 'I', currJoint.inertiaTensor);
            end
            obj.model = SerialLink(L, 'name', '4DofArm');        
        end
    end
    
    methods
        function obj = ArmModelRT(jointParams)            
            obj = obj@GeneralArmModel(jointParams);
            obj.instantiateModel();
        end
        
        function ddqVec = forwardDynamics(obj, torques, doUpdate)
            % Function has to be called only after calling update state
            state = obj.getState();
            ddqVec = obj.model.accel(state(1,:), state(2,:),...
                torques);
            
            if ~exist('doUpdate', 'var')
                doUpdate = 0;
            end
               
            if doUpdate
                obj.updateAccelerations(ddqVec);
            end
        end
        
        function torques = inverseDynamics(obj, accelerations)
            stateVec = obj.getState();
            torques = obj.model.rne(stateVec(1,:), stateVec(2,:),...
                accelerations);
        end
        
        function inertialMatrix = getInertialMatrix(obj)
            stateVec = obj.getState();
            inertialMatrix = obj.model.inerta(stateVec(1,:));
        end
        
        function x = getEndEffectorPosition(obj)
            state = obj.getState();
            trans = obj.model.fkine(state(1,:));
            [r, x] = tr2rt(trans);
        end
        
        function result = getDerivativesRelativeToState(obj, control)
            dofs = obj.totalActuatedJoints;
            result = zeros(dofs,dofs*2);
        end
    end   
 end