%%Visualizer class to make pretty pictures
classdef rlVisualizer < handle
   properties (SetAccess = private)
      models = rlCModel.empty;
      c_window_ptr = 0; 
   end
  
   
   methods
       function obj = rlVisualizer(name,width,height)
          obj.c_window_ptr = DynamicsModelMatlab(5,1,name,width,height);
       end
       
       function setBackgroundColour(obj, colour)
           if(numel(colour) ~= 3)
              error('Colour should be [R G B] array') 
           end
           DynamicsModelMatlab(5, 12, obj.c_window_ptr, colour);
       end
       
       function addModel(obj,model)
          if(~isa(model,'rlCModel') && ~find(strcmp(superclasses(model),'rlCModel')))
              error('Visualizer: Cannot add non model object'); 
          end 
          obj.models(end+1) = model;
          %Add model to the visualizer
          DynamicsModelMatlab(5,3,obj.c_window_ptr,model.c_ptr);
       end
       
       function addScene(obj,scene)
           DynamicsModelMatlab(5,8,obj.c_window_ptr,scene);
       end
       
       function add3DModel(obj,file_name,pose)
           %Adds an Assimp model to the scene, 
           %pose is defined as [x,y,z,roll,pitch,yaw]
           DynamicsModelMatlab(5,9,obj.c_window_ptr,file_name,pose);
       end
       
       function update(obj)
           DynamicsModelMatlab(5,4,obj.c_window_ptr);
       end
       
       function [I, success] = getScreenshot(obj)
           [success, I] = DynamicsModelMatlab(5, 13, obj.c_window_ptr);
           I = flipud(I);
       end
       
       function addMarker(obj, name, position, colour, visibility)
           if nargin < 4
               colour = [0; 0; 1; 0.8];
           end
           if nargin < 5
               visibility = 1;
           end
           if(numel(colour) ~= 4)
              error('Colour should be [R G B A] array') 
           end
           DynamicsModelMatlab(5, 5, obj.c_window_ptr, name, position, colour, visibility);
       end
       
       function setMarkerCovar(obj, name, cov, colour, visibility)
           % Visualize Marker Covariance
           if nargin < 4
               colour = [1; 0; 0; 0.2];
           end
           if nargin < 5
               visibility = 1;
           end
           if(numel(colour) ~= 4)
              error('Colour should be [R G B A] array') 
           end
           DynamicsModelMatlab(5, 11, obj.c_window_ptr, name, cov, colour, visibility);
       end
       
       %Visualize Sensor Covariance
       function setSensorCovar(obj,model,sensor,cov,visibility)
           DynamicsModelMatlab(5,6,obj.c_window_ptr,model.c_model_ptr,sensor.c_sensor_ptr,cov,visibility);
       end
       
       %Change Sensor Color
       function setSensorColor(obj,model,sensor,color)
           if(numel(color) ~= 4)
              error('Color should be [R G B A] array') 
           end
           DynamicsModelMatlab(5,7,obj.c_window_ptr,model.c_ptr,sensor.c_sensor_ptr,color);
       end
       
       function setTransformColor(obj,model,transform,color)
           %Change the color of a partivular transform
           
           if(numel(color) ~= 4)
               error('Color should be [R G B A] array')
           end
           if(transform.c_ptr ~= 0)
               DynamicsModelMatlab(5,16,obj.c_window_ptr,model.c_ptr,transform.c_ptr,color);
           end
       end
       
       function pose = getViewportPose(obj)
           pose = DynamicsModelMatlab(5,14,obj.c_window_ptr);
       end
       
       function setViewportPose(obj,TMat)
           DynamicsModelMatlab(5,15,obj.c_window_ptr,TMat);
       end
       
       function addCamera(obj, camera)
           DynamicsModelMatlab(5,10,obj.c_window_ptr, camera.c_ptr);
       end
       
       %Destroy c stuff
       function delete(obj)
           %delete window
           if(obj.c_window_ptr ~= 0)
                DynamicsModelMatlab(5,2,obj.c_window_ptr);
           end
       end 
   end 
end