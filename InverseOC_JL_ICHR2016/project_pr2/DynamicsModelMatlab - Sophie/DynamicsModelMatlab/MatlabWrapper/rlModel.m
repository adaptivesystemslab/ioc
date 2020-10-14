classdef rlModel < matlab.mixin.Heterogeneous & handle
    
   properties
       name = [];
       manufacturer = [];
       
       position = [];
       velocity = [];
       acceleration = [];
       torque = [];
       
       joints = rlJoint.empty;
       frames = rlFrame.empty;
       bodies = rlBody.empty;
   end
    
    
   %Protected C++ pointer
   properties (SetAccess = protected)
       %This is the actual integer pointer to C++ object in mem
       c_ptr = 0;
   end 
    
   methods
      %Constructor to load model from file
      function obj = rlModel(filename)
          %The correct loader will be called based on subclass static load
          %method, Kinemati for rlKinModel and Dynamic for rlDynModel 
          obj.c_ptr = obj.load(filename);
          
          num_bodies = DynamicsModelMatlab(1,6,obj.c_ptr);
          %Note array index in c starts at 0
          for i=0:num_bodies-1
              obj.bodies(end+1) = rlBody(DynamicsModelMatlab(1,7,obj.c_ptr,i));
              obj.bodies(end).m = 0; % mass inits to 1 for some reason
          end
          
          num_joints = DynamicsModelMatlab(1,12,obj.c_ptr);
          %Note array index in c starts at 0
          for i=0:num_joints-1
              obj.joints(end+1) = rlJoint(DynamicsModelMatlab(1,11,obj.c_ptr,i));
              switch(obj.joints(end).type)
                  case 'spherical'
                      obj.joints(end).position = [1 0 0 0];
              end
          end
      end
      
      function r = eq(obj1,obj2)
         %Returns true if the wrapper is for the same c++ object
         
         r = obj1.c_ptr == [obj2.c_ptr]; 
      end
      
      function value = get.position(obj)
          value = DynamicsModelMatlab(1,22,obj.c_ptr);
      end
      
      function set.position(obj,value)
          DynamicsModelMatlab(1,31,obj.c_ptr,value);
      end
      
      function value = get.velocity(obj)
          value = DynamicsModelMatlab(1,25,obj.c_ptr);
      end
      
      function set.velocity(obj,value)
          DynamicsModelMatlab(1,33,obj.c_ptr,value);
      end
      
      function value = get.acceleration(obj)
          value = DynamicsModelMatlab(1,5,obj.c_ptr);
      end
      
      function set.acceleration(obj,value)
          DynamicsModelMatlab(1,27,obj.c_ptr,value);
      end
      
      function value = get.torque(obj)
          value = DynamicsModelMatlab(1,23,obj.c_ptr);
      end
      
      function set.torque(obj,value)
          DynamicsModelMatlab(1,32,obj.c_ptr,value);
      end
      
   end
   
   methods(Access = protected)
      %Some protected methods used to initialize the properties on the model 
      
      
       
   end
   
    methods(Static=true,Access=protected)
      %This is a protected method which loads the model from XLM file
      %In the constructor of the object
       function c_ptr = load(filename)
            c_ptr = DynamicsModelMatlab(0);
       end
    end
    
end