classdef rlKinModel < rlModel
%Wrapper for rl::mdl::Kinematic object
    
   %Kinematic Properties of the model
   properties
       invJ = [];
       J = [];
       JdQd = [];
   end
   
   methods
      %Constructor  
      function obj = rlKinModel(filename)
          %Load the model using superclass
          obj = obj@rlModel(filename);
      end
      
      %Wrapped Methods
      
      function bool = calculateInversePosition(obj,eeTransform,eeNumber)
          %Calculate Q based on end effector transform           
          bool = DynamicsModelMatlab(2,3,obj.c_ptr,eeTransform,eeNumber);
      end
      
      function calculateJacobian(obj)
          DynamicsModelMatlab(2,4,obj.c_ptr);
      end
      
      function calculateJacobianDerivative(obj)
          DynamicsModelMatlab(2,5,obj.c_ptr);
      end
      
      function calculateJacobianInverse(obj)
          DynamicsModelMatlab(2,6,obj.c_ptr);
      end
      
      function measure = calculateManipulabilityMeasure(obj)
          measure = DynamicsModelMatlab(2,7,obj.c_ptr);
      end
      
      function forwardAcceleration(obj)
          DynamicsModelMatlab(2,8,obj.c_ptr);
      end
      
      function forwardPosition(obj)
          DynamicsModelMatlab(2,9,obj.c_ptr);
      end
      
      function forwardVelocity(obj)
          DynamicsModelMatlab(2,10,obj.c_ptr);
      end
      
      function value = get.J(obj)
          value = DynamicsModelMatlab(2,11,obj.c_ptr);
      end
      
      function value = get.JdQd(obj)
          value = DynamicsModelMatlab(2,12,obj.c_ptr);
      end
      
      function value = get.invJ(obj)
          value = DynamicsModelMatlab(2,13,obj.c_ptr);
      end
      
      %Destructor
      function delete(obj)
          if obj.c_ptr ~= 0
            DynamicsModelMatlab(2,2,obj.c_ptr);
          end
          obj.c_ptr = 0;
      end
   end 
   
   %%%%%%%%%%%%%%% SOME C++ STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   methods(Static=true,Access=protected)
      %This is a protected method which loads the model from XLM file
      %In the constructor of the object
       function c_ptr = load(filename)
            c_ptr = DynamicsModelMatlab(2,1,filename);
       end
    end
end