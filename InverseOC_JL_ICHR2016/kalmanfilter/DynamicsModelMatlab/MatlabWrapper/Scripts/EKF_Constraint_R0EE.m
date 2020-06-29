classdef EKF_Constraint_T0EE < EKF_Constraint
   %EKF Constraint such that the transformation from end effector to base frame
   %matches. For example consider holding a bar in two hands, for both end
   %effectors the bar goes along Y direction in end effector frame. Thus
   %the projection of Y of each end effector onto the base frame should be
   %the same. Thus the Y part (second column) of R0EE should be the same.
   
   properties       
       %Type selects the axes that should match, default is all enabled for
       %the bar example described above type would be [0 1 0]
       type = true(3,1);
   end
   
   methods
       
       function obj = EKF_Constraint_R0EE(name,model1,frame1,T1,model2,frame2,T2)
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
          dR0EE = obj.sensors(1).dTb_ee(1:3,1:3,:)
           J1 = obj.sensors(1).jacobian(obj.type,:);
           J2 = obj.sensors(2).jacobian(obj.type,:);
          value = [J1, -J2];
       end
       
       function value = get_b_value(obj)
           %b part of Ax=b constraint
           
           T1 = obj.sensors(1).transform;
           T2 = obj.sensors(2).transform;
           [yd, pd, rd] = dcm2angle((T2(1:3,1:3)*T1(1:3,1:3)')');
           rd = [yd;pd;rd];
           t1 = T1(1:3,4);
           t2 = T2(1:3,4);
           value = [rd;t2-t1] + ...
               obj.A * [obj.models(1).position; obj.models(2).position];
           value = value(obj.type);
       end
   end
   
    
end