classdef rlJoint < rlTransform
      
    properties
        position;
        velocity;
        acceleration;
    end
    
   methods 
       
    function obj = rlJoint(c_ptr)
          %C++ pointer stuff
          obj = obj@rlTransform(c_ptr);
          %Get the type of joint
          obj.type = DynamicsModelMatlab(10,1,obj.c_ptr);
    end
    
    function value = get.position(obj)
       value = DynamicsModelMatlab(10,2,obj.c_ptr); 
    end
    
    function set.position(obj,value)
        DynamicsModelMatlab(10,3,obj.c_ptr,value); 
    end
    
    function value = get.velocity(obj)
       value = DynamicsModelMatlab(10,4,obj.c_ptr); 
    end
    
   function set.velocity(obj,value)
        DynamicsModelMatlab(10,5,obj.c_ptr,value); 
    end
    
    function value = get.acceleration(obj)
       value = DynamicsModelMatlab(10,6,obj.c_ptr); 
    end
    
    function set.acceleration(obj,value)
        DynamicsModelMatlab(10,7,obj.c_ptr,value); 
    end
    
   end
    
    
end