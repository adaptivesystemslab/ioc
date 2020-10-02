classdef EKF_Resizable < EKF
   %Here we have EKF based on a kinematic model that is automatical resized
   %if sensors are missing 
   properties(SetAccess = protected)
      %Flag if the sensor is missing or not (true = visible, false = missing)
      visible_indexes = logical.empty();
      visible_sensors = logical.empty();
      full_observation = [];
      model_handle = [];  
   end
   
   methods
      function obj = EKF_Resizable(model)
           obj.model_handle = model; 
           
           obj.visible  = true(numel(model.sensors),1);
           %Initialize Stuff
           z = vertcat(obj.model_handle.sensors.measurement);
           obj.observation_noise = eye(numel(z));
           obj.process_noise = eye(model.dof*2);
           obj.covariance = eye(model.dof*2);
           obj.state = [model.position;model.velocity];
           
           obj.inov_covariance = zeros(numel(z),numel(z)); 
           
           obj.sizeZ = numel(z);
           obj.sizeX = size(obj.process_noise,1);
      end
      
      
        
      function z = makeVisMeasure(obj,x)
            %Measurement is based on the sensors attached to the model  
            z = obj.makeMeasure(x);
            z = z(obj.visible_indexes);
      end
      
      function H = makeH(obj,x)
          %Make the jacobian for visible sensors, can be overridden 
          H = vertcat(obj.model_handle.sensors(obj.visible).obsJacobian);
      end

      function H = makeFullH(obj,x)
          %Make the jacobian for all sensors, can be overridden 
         H = vertcat(obj.model_handle.sensors.obsJacobian); 
      end
      
      function H = makeMissingH(obj,x)
          %Make the jacobian for missing markers, can be overridden
         H = vertcat(obj.model_handle.sensors(~obj.visible).obsJacobian); 
         if isempty(H)
            H = zeros(0,obj.model_handle.dof*3); 
         end
      end
      
      function probs = findMatches(obj,measurements)
         %Computes the probability of each measurement to belong to one of the 
         %missing sensors. Measurements should be an array where each
         %row represents a potential measurement for the sensors
         
         %Get the missing marker jacobian
         H = obj.makeMissingH(obj.state);
         
         %Project covariance 
         innov_all = H*obj.covariance*H';
         
         missing_sensors = obj.model_handle.sensors(~obj.visible);
         
         %Init 
         probs = zeros(numel(missing_sensors),size(measurements,1));
         
         indx = 1;
         
         for i=1:numel(missing_sensors)
             sens = missing_sensors(i);
             %Make row vector
             sens_mes = sens.measurement';
             %Grab idependend error covariance
             sigma = innov_all(indx:indx+numel(sens_mes)-1,indx:indx+numel(sens_mes)-1);
             
             %Grab only part of covariance
             if size(measurements,2) < numel(sens_mes)
                sigma = sigma(1:size(measurements,2),1:size(measurements,2));
                sens_mes = sens_mes(1:size(measurements,2));
             elseif size(measurements,2) > numel(sens_mes)
                 error('Potential measurement should not be a larger vector than sensor measurement.');
             end
             %trololololo hacky psd
             sigma = ( 1/2.*sigma + 1/2.*sigma' );
             probs(i,:) = mvnpdf(measurements,sens_mes,sigma)./mvnpdf(sens_mes,sens_mes,sigma);
         end
      end
      
      function x_new = stateUpdate(obj,x,u)
           %State Update, integrate Velocity into positions where u is dt
           dof = obj.model_handle.dof;
           x_new = x;
           x_new(1:dof) = x_new(1:dof)+...
               x(dof+1:end)*u;
           
            %Forward Kinematics
            dof = obj.model_handle.dof;
            obj.model_handle.position = x_new(1:dof);
            obj.model_handle.velocity = x_new(dof+1:end);
            obj.model_handle.forwardKinematics();
           
        end
        
        function A = makeA(obj,x,u)
            dof = obj.model_handle.dof;
            A = eye(numel(x),numel(x)); 
            A(1:dof,dof+1:end) = A(1:dof,dof+1:end)+eye(dof,dof)*u;
            A(1:dof,dof*2+1:end) = A(1:dof,dof*2+1:end)+eye(dof,dof)*u*u/2;
        end
   end
   
   
   methods (Access=protected)
       
       function value = MakeObservationNoise(obj)
            %Only get observation noise from not missing sensors
            mes_sizes = cellfun(@numel,{obj.model_handle.sensors.measurement});
            indexes = cumsum(mes_sizes) - mes_sizes(1)+1;
            
            old_noise = diag(obj.observation_noise);
            noise = zeros(sum(mes_sizes(obj.visible)),1);
            indx_small = 1;
            for i=1:numel(mes_sizes)
                if obj.visible(i)
                    noise(indx_small:indx_small+mes_sizes(i)-1) = old_noise(indexes(i):indexes(i)+mes_sizes(i)-1);
                    indx_small = indx_small + mes_sizes(i);
                end
            end
            
            value = diag(noise);
            
       end
       
   end
    
end