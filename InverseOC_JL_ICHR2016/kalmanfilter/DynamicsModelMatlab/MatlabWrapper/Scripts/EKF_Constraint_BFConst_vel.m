classdef EKF_Constraint_BFConst_vel <  EKF_Constraint_BF_vel
    properties
       %Desired velocity in world frame [w v]';
       Dv = zeros(6,1);
    end
    
    properties(SetAccess=protected)
        trans;
    end
    
    
    methods
       function obj = EKF_Constraint_BFConst_vel(name,model,frame,T,Dv)
          %Constraint wrt constant transform
          obj = obj@EKF_Constraint_BF_vel(name,model,frame,T,model,'world',Dv); 
          obj.Dv = Dv;
          %Pull out the transform of the second sensor (one that is attached to world)
          t_indx = find(contains({model.transforms.name},[name '_Bt']),1);
          obj.trans = model.transforms(t_indx);
       end 
       
       function set.Dv(obj,value)
           %Set the constraint in world frame 
           obj.Dv = value;
       end
       
    end
    
end