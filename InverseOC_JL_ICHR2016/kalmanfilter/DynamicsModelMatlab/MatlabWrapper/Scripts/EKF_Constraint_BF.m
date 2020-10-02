classdef EKF_Constraint_BF < EKF_Constraint
    %EKF constraint in base frame
    properties
        %Type selects the rows of base frame Jacobians used. For example if
        %we want to only constrain rotation about Z and position in Y we set
        %type to [0 0 1 0 1 0]
        type = true(6,1);
    end
    
    methods
        
        function obj = EKF_Constraint_BF(name,model1,frame1,T1,model2,frame2,T2)
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
            
            T1 = obj.sensors(1).transform;
            T2 = obj.sensors(2).transform;
            [yd, pd, rd] = dcm2angle2((T2(1:3,1:3)*T1(1:3,1:3)')'); % dcm2angle is in aerospace toolbox. copied out the function so we don't invoke it, UW doesn't have a lot of license for this
            rd = [rd;pd;yd];
            t1 = T1(1:3,4);
            t2 = T2(1:3,4);
            value = [rd;t2-t1];
            if obj.models(1) == obj.models(2)
                value = value(obj.type) + obj.A * obj.models(1).position;
            else
                value = value(obj.type) + ...
                    obj.A * [obj.models(1).position; obj.models(2).position];
            end
        end
    end
    
    
end