classdef rlCModel < rlDynModel
   %This is the extension of robotics library that allows Jacobian
   %calculation and stuff
   
   properties
      sensors = SensorCore.empty;
      base = [];
      transforms = rlTransform.empty;
      
      %Magnetic Field
      m;
   end
   
   methods
      
       function obj = rlCModel(filename)
           %Load the model using superclass
           obj = obj@rlDynModel(filename);
           
           num_transforms = DynamicsModelMatlab(4,11,obj.c_ptr);
           
           joint_ptrs = [obj.joints.c_ptr];
           
           for i=0:num_transforms-1
               %We dont want to add joints to our transform vector even
               %through they are transfroms or do we ?
               c_ptr = DynamicsModelMatlab(4,12,obj.c_ptr,i);
               if(isempty(find(joint_ptrs == c_ptr,1)))
                    obj.transforms(end+1) = rlTransform(c_ptr);
               else
                    obj.transforms(end+1) = rlJoint(c_ptr);
               end
           end
       end
       
       function addSensor(obj,sensor,frame,transform)
           
           existing_indx = find(contains({obj.sensors.name},sensor.name),1);
           if ~isempty(existing_indx)
              warning(['Model ' obj.name ' already contains sensor ' sensor.name ' it will be removed and added again']);
              obj.removeSensor(obj.sensors(existing_indx));
           end
           
           
            DynamicsModelMatlab(4,4,obj.c_ptr,sensor.c_sensor_ptr,frame,transform);
            obj.sensors(end+1) = sensor;
            %Get the transform into our transforms vector, we have to loop
            %through all of the transforms 
            num_transforms = DynamicsModelMatlab(4,11,obj.c_ptr);
            obj.transforms = rlTransform.empty;
            joint_ptrs = [obj.joints.c_ptr];
            for i=0:num_transforms-1
                %We dont want to add joints to our transform vector even
                %through they are transfroms or do we ?
                c_ptr = DynamicsModelMatlab(4,12,obj.c_ptr,i);
                if(isempty(find(joint_ptrs == c_ptr,1)))
                    obj.transforms(end+1) = rlTransform(c_ptr);
                else
                    obj.transforms(end+1) = rlJoint(c_ptr);
                end
            end
       end
       
       function bool = removeSensor(obj,sensor)
           bool = false;
           sensor_indx = find(sensor == obj.sensors);
           if ~isempty(sensor_indx)
              obj.sensors(sensor_indx) = [];
              bool = DynamicsModelMatlab(4,5,obj.c_ptr,sensor.c_sensor_ptr);
           end
       end
        
       function calculateSensorJacobians(obj)
           DynamicsModelMatlab(4,14,obj.c_ptr);
       end
       
       function set.m(obj,value)
            DynamicsModelMatlab(4,19,obj.c_ptr,value); 
       end
       
      function value = get.m(obj)
            value = DynamicsModelMatlab(4,16,obj.c_ptr); 
       end
        
        function set.base(obj,frame_name)
           DynamicsModelMatlab(4,9,obj.c_ptr,frame_name); 
           obj.base = frame_name;
           
           num_transforms = DynamicsModelMatlab(4,11,obj.c_ptr);
           joint_ptrs = [obj.joints.c_ptr];
           obj.transforms = rlTransform.empty;
           for i=0:num_transforms-1
               %We dont want to add joints to our transform vector even
               %through they are transfroms or do we ?
               c_ptr = DynamicsModelMatlab(4,12,obj.c_ptr,i);
               if(isempty(find(joint_ptrs == c_ptr,1)))
                    obj.transforms(end+1) = rlTransform(c_ptr);
               else
                    obj.transforms(end+1) = rlJoint(c_ptr);
               end
           end
        end
        
        function frame = getFrameByName(obj,name)
            frame = [];
            frame_cptr = DynamicsModelMatlab(4,13,obj.c_ptr,name);
            if frame_cptr ~= 0
                frame = rlFrame(frame_cptr);
            end
        end
        
       %Destructor
      function delete(obj)
          %Destructor
          if obj.c_ptr ~= 0
            DynamicsModelMatlab(4,2,obj.c_ptr);
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
            c_ptr = DynamicsModelMatlab(4,1,filename);
       end
    end
    
end