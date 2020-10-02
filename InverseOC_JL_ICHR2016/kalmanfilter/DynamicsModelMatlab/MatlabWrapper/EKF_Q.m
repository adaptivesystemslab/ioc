classdef EKF_Q < EKF_Resize
    
    properties
         
    end
    
    methods
        
        function obj = EKF_Q(model)
            
            %Call super
            obj = obj@EKF_Resize(model);
            
            %Our observations are the sensors attached to the model
            z = vertcat(model.sensors.measurement);
            obj.observation_noise = eye(numel(z));
            %Our state is model q dq ddq
            obj.state = [model.position];
            obj.process_noise = eye(numel(obj.state));
            obj.covariance = diag(numel(obj.state));
            obj.inov_covariance = zeros(numel(z),numel(z));            
            obj.sizeZ = numel(z);
            obj.sizeX = numel(obj.state);
        end
        
          
        function z = makeMeasure(obj,x)
            %Predict measurement from state
            dof = numel(obj.model_handle.joints);
            obj.model_handle.position = x(1:dof);
            obj.model_handle.forwardPosition();
            %Generate measurements
            z = SensorMeasurement(obj.model_handle.sensors);

        end
        
        function x_new = stateUpdate(obj,x,u)
               %State Update
               
               dof = numel(obj.model_handle.joints);
               
               %Constant Q
               x_new = x;
                
               %Forward Kinematics
               obj.model_handle.position = x_new(1:dof);
               obj.model_handle.forwardPosition();
        end
        
        function A = makeA(obj,x,u)
            A = eye(numel(x),numel(x));
        end
        
       function H = makeH(obj,x)
            %Generate Sensor Jacobian Matrix based on sensors attached to
            %the model
            obj.model_handle.calculateSensorJacobians();
            H = vertcat(obj.model_handle.sensors.obsJacobian);
            H = H(:,1:numel(x));
       end
    end        
end