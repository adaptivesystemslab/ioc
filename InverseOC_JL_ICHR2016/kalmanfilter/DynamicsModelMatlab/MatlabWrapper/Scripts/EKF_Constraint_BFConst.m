classdef EKF_Constraint_BFConst <  EKF_Constraint_BF
    properties
       Tw = eye(4); 
    end
    
    properties(SetAccess=protected)
        trans;
    end
    
    
    methods
       function obj = EKF_Constraint_BFConst(name,model,frame,T,Tw)
          %Constraint wrt constant transform
          obj = obj@EKF_Constraint_BF(name,model,frame,T,model,'world',Tw); 
          obj.Tw = Tw;
          %Pull out the transform of the second sensor (one that is attached to world)
          t_indx = find(contains({model.transforms.name},[name '_Bt']),1);
          obj.trans = model.transforms(t_indx);
       end 
       
       function set.Tw(obj,value)
           %Set the constraint in world frame 
           obj.Tw = value;
           obj.trans.t = obj.Tw;
       end
       
    end
    
end