classdef Body < Frame
   properties
       fX = zeros(6,1);
       com = zeros(3,1);
       I = eye(3,3);
       m = 0;
   end
   
   methods
       
       function obj = Body(c_body_ptr)
          obj = obj@Frame(c_body_ptr); 
       end
       
       function set.com(obj,value)
           DynamicsModelMatlab(8,1,obj.c_frame_ptr,value); 
       end
       
       function set.m(obj,value)
           DynamicsModelMatlab(8,2,obj.c_frame_ptr,value); 
       end
       
       function set.I(obj,value)
           if size(value,1) == 3 && size(value,2) == 3
               %Set inertia requires [xx, yy, zz, yz, xz, xy]
                DynamicsModelMatlab(8,3,obj.c_frame_ptr,[value(1,1) value(2,2) value(3,3) value(2,3) value(1,3) value(1,2)]);
           else
                DynamicsModelMatlab(8,3,obj.c_frame_ptr,value);
           end
       end
       
       function value = get.com(obj)
          value  = DynamicsModelMatlab(8,4,obj.c_frame_ptr); 
       end
       
       function value = get.fX(obj)
          value  = DynamicsModelMatlab(8,5,obj.c_frame_ptr); 
       end
       
       function value = get.I(obj)
          value  = DynamicsModelMatlab(8,6,obj.c_frame_ptr); 
       end
       
       function value = get.m(obj)
          value  = DynamicsModelMatlab(8,7,obj.c_frame_ptr); 
       end
       
       function set.fX(obj,value)
           DynamicsModelMatlab(8,8,obj.c_frame_ptr,value);
       end
   end 
end