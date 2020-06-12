classdef SensorCore < handle
   %This is the Core Sensor which actually creates the body which gets
   %attached to the Model.
    properties(SetAccess=protected)
       c_sensor_ptr = [];
   end
   
   properties (SetAccess=protected)
       velocity = [0 0 0 0 0 0];
       acceleration = [0 0 0 0 0 0];
       %Observation Jacobian which changes based on decorators
       obsJacobian = [];
       %Regular Jacobian (end effector frame)
       jacobian = [];
       %Base Frame Jacobian
       baseJacobian = [];
       %Base Frame Jacobian Derivative
       baseJacobianDerivative = [];
       %Derivative of base to ee transform wrt each joint
       dTb_ee = [];
       %Derivative of ee to base transform wrt each joint
       dTee_b = [];
       %Derivative of angular velocity in sensor frame wrt joints
       dWdQ = [];
       
       measurement = [];
       transform = [];
       name = [];
       type = [];
       binary_type = 0;
   end
   
   properties
      %The base frame of the sensor, world by default
       base = []; 
   end
   
   methods
       
       function obj = SensorCore(name)
           %Create the core sensor which is the actual Body which attached to the Model
           obj.c_sensor_ptr = DynamicsModelMatlab(6,1,name);
           obj.name = name;
        end
        
        function addDecorator(obj,type)
           %Add additional decorator to sensor
           switch type
               case 'gyroscope'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'gyroscope');
               case 'accelerometer'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'accelerometer');
               case 'magnetometer'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'magnetometer');
               case 'orientation'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'orientation');
               case 'position'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'position');
               case 'velocity'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'velocity');
               case 'angularvelocity'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'angularvelocity');
               case 'yaw'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'yaw');
               case 'constbfpx'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'constbfpx');
               case 'constbfpy'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'constbfpy');
               case 'constbfpz'
                   obj.c_sensor_ptr = DynamicsModelMatlab(6,3,obj.c_sensor_ptr,'constbfpz');
               otherwise
                   warning('Unknown SensorDecorator Type -> IGNORING');
           end
        end
        
        function value = get.type(obj)
           value =  DynamicsModelMatlab(6,9,obj.c_sensor_ptr);
        end
        
        function value = get.binary_type(obj)
           value =  DynamicsModelMatlab(6,18,obj.c_sensor_ptr);
        end
           
       function value = get.velocity(obj)
       % Returns the linear and angular velocity of the sensor in 
       % Sensor Frame
           obj.velocity = DynamicsModelMatlab(6,4,obj.c_sensor_ptr);
           value = obj.velocity;
       end
       
       function value = get.acceleration(obj)
       % Returns the linear and angular acceleration of the sensor in 
       % Sensor Frame
           obj.acceleration = DynamicsModelMatlab(6,5,obj.c_sensor_ptr);
           value = obj.acceleration;
       end
       
       
       function value = get.measurement(obj)
       % Returns the linear and angular acceleration of the sensor in 
       % Sensor Frame
           obj.measurement = DynamicsModelMatlab(6,6,obj.c_sensor_ptr);
           value = obj.measurement;
       end
       
       function value = get.obsJacobian(obj)
           obj.obsJacobian = DynamicsModelMatlab(6,7,obj.c_sensor_ptr);
           value = obj.obsJacobian;
       end
       
       function value = get.jacobian(obj)
           obj.jacobian = DynamicsModelMatlab(6,11,obj.c_sensor_ptr);
           value = obj.jacobian;
       end
       
       function value = get.baseJacobian(obj)
           obj.baseJacobian = DynamicsModelMatlab(6,12,obj.c_sensor_ptr);
           value = obj.baseJacobian;
       end
       
       function value = get.baseJacobianDerivative(obj)
           obj.baseJacobianDerivative = DynamicsModelMatlab(6,19,obj.c_sensor_ptr);
           value = obj.baseJacobianDerivative;
       end
       
      function value = get.transform(obj)
           obj.transform = DynamicsModelMatlab(6,8,obj.c_sensor_ptr);
           value = obj.transform;
      end
       
      function value = get.dTb_ee(obj)
          value = DynamicsModelMatlab(6,13,obj.c_sensor_ptr);
      end
      
      function value = get.dTee_b(obj)
          value = DynamicsModelMatlab(6,14,obj.c_sensor_ptr);
      end
      
      function value = get.dWdQ(obj)
         value =  DynamicsModelMatlab(6,17,obj.c_sensor_ptr);
      end
      
      function set.base(obj,frame)
         if ~isa(frame,'rlFrame')
            error('Base of sensor has to be a valid frame object'); 
         end
         DynamicsModelMatlab(6,20,obj.c_sensor_ptr,frame.c_frame_ptr);
      end
      
      function value = get.base(obj)
          c_frame_ptr = DynamicsModelMatlab(6,21,obj.c_sensor_ptr);
          if c_frame_ptr ~= 0
             value = rlFrame(c_frame_ptr);
          else
              value = [];
          end
      end
      
      function iseq = eq(obj1,obj2)
         iseq =  ([obj1.c_sensor_ptr] == [obj2.c_sensor_ptr])';
      end
       
       function delete(obj)
           %delete c sensor
            if(obj.c_sensor_ptr ~= 0)
                %Sensor C pointer is deleted when the model is deleted
                %This cleans up all the decorators
                DynamicsModelMatlab(6,2,obj.c_sensor_ptr);
            end 
       end
   end
end