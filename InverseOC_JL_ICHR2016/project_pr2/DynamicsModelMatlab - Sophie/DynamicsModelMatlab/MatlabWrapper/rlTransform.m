classdef rlTransform < handle & matlab.mixin.Heterogeneous 
   
    properties
       t = []; 
    end
    
    properties (SetAccess = protected)
       %Protected properties
       name = '';
       type = 'fixed';
       c_ptr = 0;
       frame_in = [];
       frame_out = [];
       x = [];
    end
    

    
    methods
        function obj = rlTransform(c_ptr)
            obj.c_ptr = c_ptr;
        end
        
        function value = get.name(obj)
            value = DynamicsModelMatlab(9,5,obj.c_ptr);
        end
        
        function value = get.frame_in(obj)
           value = rlFrame(DynamicsModelMatlab(9,1,obj.c_ptr)); 
        end
        
        function value = get.frame_out(obj)
           value = rlFrame(DynamicsModelMatlab(9,2,obj.c_ptr)); 
        end
        
        function value = get.t(obj)
           value = DynamicsModelMatlab(9,3,obj.c_ptr); 
        end
        
        function set.t(obj,value)
           DynamicsModelMatlab(9,6,obj.c_ptr,value); 
        end
        
        function value = get.x(obj)
           value  = DynamicsModelMatlab(9,4,obj.c_ptr);
       end
    end
end