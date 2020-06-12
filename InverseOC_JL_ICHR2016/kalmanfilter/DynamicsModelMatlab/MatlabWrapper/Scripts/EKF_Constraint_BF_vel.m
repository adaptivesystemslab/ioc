classdef EKF_Constraint_BF_vel < EKF_Constraint
    %EKF constraint in base frame
    properties
        %Type selects the rows of base frame Jacobians used. For example if
        %we want to only constrain rotation about Z and position in Y we set
        %type to [0 0 1 0 1 0]
        type = true(6,1);
    end
    
    methods
        
        function obj = EKF_Constraint_BF_vel(name,model1,frame1,T1,model2,frame2,T2)
            %Constructor where all the tings are defined
            
            obj = obj@EKF_Constraint(name,model1,frame1,T1,model2,frame2,T2);
        end
    end
    
    methods(Access=protected)
        function value = get_A_value(obj)
            %A part of Ax=b constraint, this is the sensor Jacobians
            
            %Run calculate jacobians if needed
            if(isempty(obj.sensors(1).baseJacobian))
                obj.models(1).calculateSensorJacobians();
            end
            if(isempty(obj.sensors(2).baseJacobian))
                obj.models(2).calculateSensorJacobians();
            end
            
            %Build A
            J1 = obj.sensors(1).baseJacobian(obj.type,:);
            J2 = obj.sensors(2).baseJacobian(obj.type,:);
            if obj.models(1) == obj.models(2)
                value = [J1-J2];
            else
                value = [J1, -J2];
            end
        end
        
        function value = get_b_value(obj)
            %b part of Ax=b constraint
            R = obj.sensors(1).transform(1:3,1:3);
            R = blkdiag(R,R);
            V1 = R*obj.sensors(1).velocity;
            V2 = R*obj.sensors(2).velocity;
            value = V2 - V1;
            if obj.models(1) == obj.models(2)
                value = value(obj.type) + obj.A * obj.models(1).velocity;
            else
                value = value(obj.type) + ...
                    obj.A * [obj.models(1).velocity; obj.models(2).velocity];
            end
        end
    end
    
    
end