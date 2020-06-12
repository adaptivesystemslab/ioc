classdef EKF_Constraint_T0EE < EKF_Constraint
    %EKF Constraint such that the transformation from end effector to base frame
    %matches. For example consider holding a bar in two hands, for both end
    %effectors the bar goes along Y direction in end effector frame. Thus
    %the projection of Y of each end effector onto the base frame should be
    %the same. Thus the Y part (second column) of R0EE should be the same.
    %please note that when position constraints are enabled you should also
    %enable the rotation constraints.
    
    properties
        %Type selects the axes that should match, default is all enabled for
        %the bar example described above type would be [0 1 0 0 0 0].
        type = true(6,1);
    end
    
    methods
        
        function obj = EKF_Constraint_T0EE(name,model1,frame1,T1,model2,frame2,T2)
            %Constructor where all the tings are defined
            obj = obj@EKF_Constraint(name,model1,frame1,T1,model2,frame2,T2);
        end
    end
    
    methods(Access=protected)
        
        function value = get_A_rot(obj)
            %Build A_rot
            
            T01 = obj.sensors(1).transform;
            T02 = obj.sensors(2).transform;
            T10 = SE3.fastInverse(T01);
            %T20 = SE3.fastInverse(T02);
            T12 = T10*T02;
            
            dT01 = obj.sensors(1).dTb_ee;
            dT10 = obj.sensors(1).dTee_b;
            dT02 = obj.sensors(2).dTb_ee;
            
            %First build the rotation constraints
            A1_rot = permute(dT01(1:3,obj.type(1:3),:),[1,3,2]);
            A1_rot = reshape(permute(A1_rot,[2,1,3]),size(A1_rot,2),[])';
            A2_rot = permute(dT02(1:3,obj.type(1:3),:),[1,3,2]);
            A2_rot = reshape(permute(A2_rot,[2,1,3]),size(A2_rot,2),[])';
            value = [A1_rot, -A2_rot];
            %Collapse if it is the same model
            if obj.models(1) == obj.models(2)
                value = A1_rot-A2_rot;
            end
        end
        
        function value = get_A_pos(obj)
            %Build A_rot
            T01 = obj.sensors(1).transform;
            T02 = obj.sensors(2).transform;
            T10 = SE3.fastInverse(T01);
            %T20 = SE3.fastInverse(T02);
            T12 = T10*T02;
            
            dT01 = obj.sensors(1).dTb_ee;
            dT10 = obj.sensors(1).dTee_b;
            dT02 = obj.sensors(2).dTb_ee;
            
            dT10_mul_T02 = reshape(permute(dT10,[2 1 3]),size(dT10,2),[])'*T02;
            %dT10_mul_T02 = permute(reshape(dT10_mul_T02',4,4,[]),[2 1 3]);
            T10_mul_dT02 = (T10*reshape(permute(dT02,[2 1 3]),size(dT02,2),[]))';
            %T10_mul_dT02 = permute(reshape(T10_mul_dT02',4,4,[]),[2 1 3]);
            %Append zeros if the constraint is between two different models
            if(obj.models(1) ~= obj.models(2))
                tmp = size(dT10_mul_T02,3);
                dT10_mul_T02 = [dT10_mul_T02; zeros(4*numel(obj.models(2).joints),4)];
                %dT10_mul_T02 = cat(3,dT10_mul_T02,zeros(4,4,size(T10_mul_dT02,3)));
                T10_mul_dT02 = [zeros(4*numel(obj.models(1).joints),4); T10_mul_dT02];
                %T10_mul_dT02 = cat(3,zeros(4,4,tmp),T10_mul_dT02);
            end
            
            dTee1_ee2 = permute(reshape((dT10_mul_T02 + T10_mul_dT02)',4,4,[]),[2 1 3]);
            dtee1_ee2 = permute(dTee1_ee2(1:3,4,:),[1 3 2]);
            value = dtee1_ee2(obj.type(4:6),:);
        end
        
        function value = get_A_value(obj)
            %A part of Ax=b constraint, this is the sensor Jacobians
            
            %Run calculate jacobians if needed
            if(isempty(obj.sensors(1).baseJacobian))
                obj.models(1).calculateSensorJacobians();
            end
            if(isempty(obj.sensors(2).baseJacobian))
                obj.models(2).calculateSensorJacobians();
            end
            
            value = [obj.get_A_rot; obj.get_A_pos];
            
        end
        
        function value = get_b_value(obj)
            %b part of Ax=b constraint
            T1 = obj.sensors(1).transform;
            T2 = obj.sensors(2).transform;
            t12 = SE3.fastInverse(T1)*T2*[0 0 0 1]';
            t12 = t12(1:3);
            R01 = T1(1:3,1:3);
            R02 = T2(1:3,1:3);
            
            if obj.models(1) ~= obj.models(2)
                if ~isempty(find(obj.type(1:3),1))
                    ax_diff = R02(:,obj.type(1:3))-R01(:,obj.type(1:3)); 
                    b_rot = ax_diff(:) +...
                        obj.get_A_rot * [obj.models(1).position; obj.models(2).position];
                else
                    b_rot = [];
                end
                b_pos = -t12(obj.type(4:6)) + obj.get_A_pos * [obj.models(1).position; obj.models(2).position];
            else
                if ~isempty(find(obj.type(1:3),1))
                    b_rot = R02(:,obj.type(1:3))-R01(:,obj.type(1:3)) +...
                        obj.get_A_rot * [obj.models(1).position];
                else
                   b_rot = []; 
                end
                b_pos = -t12(obj.type(4:6)) + obj.get_A_pos*obj.models(1).position;
            end
            value = [b_rot;b_pos];
        end
    end
    
    
end