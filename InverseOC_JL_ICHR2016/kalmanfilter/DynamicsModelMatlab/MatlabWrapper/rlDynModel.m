classdef rlDynModel < rlKinModel
    
    properties
        %Centrifugal Coriolis Vector
        V = [];
        %Gravity Vector G(q)
        G = [];
        %Mass Matrix Inverse
        invM = [];
        %Mass Matrix
        M = [];
        %Operational Mass Matrix Inverse
        invMx = [];
        %World gravity, usually [0 0 9.81]
        g
    end
    
    methods
      
      function obj = rlDynModel(filename)
          %Load the model using superclass
          obj = obj@rlKinModel(filename);
      end
      
      function calculateCentrifugalCoriolis(obj)
          DynamicsModelMatlab(3,3,obj.c_ptr);
      end
      
      function calculateGravity(obj)
          DynamicsModelMatlab(3,4,obj.c_ptr);
      end
        
      function calculateMassMatrix(obj)
          DynamicsModelMatlab(3,5,obj.c_ptr);
      end 
      
      function calculateMassMatrixInverse(obj)
          DynamicsModelMatlab(3,6,obj.c_ptr);
      end
      
      function calculateOperationalMassMatrixInverse(obj)
          DynamicsModelMatlab(3,7,obj.c_ptr);
      end
      
      function eulerCauchy(obj,dt)
          DynamicsModelMatlab(3,8,obj.c_ptr,dt);
      end
      
      function forwardDynamics(obj)
          DynamicsModelMatlab(3,9,obj.c_ptr);
      end
      
      function value = get.V(obj)
         value =  DynamicsModelMatlab(3,10,obj.c_ptr);
      end
      
      function value = get.G(obj)
          value =  DynamicsModelMatlab(3,11,obj.c_ptr);
      end
      
      function value = get.invM(obj)
          value =  DynamicsModelMatlab(3,12,obj.c_ptr);
      end
      
      function value = get.M(obj)
          value =  DynamicsModelMatlab(3,13,obj.c_ptr);
      end
      
      function value = get.invMx(obj)
          value =  DynamicsModelMatlab(3,14,obj.c_ptr);
      end
      
      function value = get.g(obj)
          value =  DynamicsModelMatlab(3,15,obj.c_ptr);
      end
      
      function inverseDynamics(obj)
          DynamicsModelMatlab(3,16,obj.c_ptr);
      end
      
      function inverseForce(obj)
           DynamicsModelMatlab(3,17,obj.c_ptr);
      end
      
      function rungeKuttaNystrom(obj,dt)
          DynamicsModelMatlab(3,18,obj.c_ptr,dt);
      end
      
      function set.g(obj,value)
          DynamicsModelMatlab(3,19,obj.c_ptr,value);
      end
      
      %Destructor
      function delete(obj)
          if obj.c_ptr ~= 0
            DynamicsModelMatlab(3,2,obj.c_ptr);
          end
          %This is needed so superclasses know not to call C++ destroctors
          %on already deleted objects
          obj.c_ptr = 0;
      end
        
    end
    
   %%%%%%%%%%%%%%% SOME C++ STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   methods(Static=true,Access=protected)
      %This is a protected method which loads the model from XLM file
      %In the constructor of the object
       function c_ptr = load(filename)
            c_ptr = DynamicsModelMatlab(3,1,filename);
       end
    end

end