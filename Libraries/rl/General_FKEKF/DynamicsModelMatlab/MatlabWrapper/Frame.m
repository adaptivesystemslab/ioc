classdef Frame < handle
   properties(SetAccess = protected)
       
       %These are the properties of a frame as defined in robotics library
       
       %Motion Vector Acceleration
       a = [];
       %Motion Vector C
       c = [];
       %Force acting on the frame
       f = [];
       %Rigid Body Inertia
       i = [];
       %Articulated Body Inertia
       iA = [];
       %Force Vector pA
       pA = [];
       %Transform of the frame
       t = [];
       %Velocity of the frame
       v = [];
       %Pleuker Transform of the frame
       x = [];
   end
   
   properties (SetAccess = protected)
       %Protected properties
       name = '';
       c_frame_ptr = 0;
   end
   
   methods 
       function obj = Frame(c_frame_ptr)
           obj.c_frame_ptr = c_frame_ptr;
       end
       
       function value = get.name(obj)
           value  = DynamicsModelMatlab(7,1,obj.c_frame_ptr);
       end
       
       function value = get.a(obj)
           value   = DynamicsModelMatlab(7,2,obj.c_frame_ptr);
       end
       
       function value = get.c(obj)
           value  = DynamicsModelMatlab(7,3,obj.c_frame_ptr);
       end
       
       function value = get.f(obj)
           value  = DynamicsModelMatlab(7,4,obj.c_frame_ptr);
       end
       
       function value = get.i(obj)
           value  = DynamicsModelMatlab(7,5,obj.c_frame_ptr);
       end
       
       function value = get.iA(obj)
           value  = DynamicsModelMatlab(7,6,obj.c_frame_ptr);
       end
       
       function value = get.pA(obj)
           value  = DynamicsModelMatlab(7,7,obj.c_frame_ptr);
       end
       
       function value = get.t(obj)
           value  = DynamicsModelMatlab(7,8,obj.c_frame_ptr);
       end
       
       function value = get.v(obj)
           value  = DynamicsModelMatlab(7,9,obj.c_frame_ptr);
       end
       
       function value = get.x(obj)
           value  = DynamicsModelMatlab(7,10,obj.c_frame_ptr);
       end
   end  
end