classdef EKF_Q_DQ < EKF_Resize
    
    properties
        
    end
    
    methods
        
        function obj = EKF_Q_DQ(model)
            
            %Call super
            obj = obj@EKF_Resize(model);
            
            %Our observations are the sensors attached to the model
            z = vertcat(model.sensors.measurement);
            obj.observation_noise = eye(numel(z));
            %Our state is model q dq ddq
            obj.state = [...
                model.position ; ...
                model.velocity ; ...
                ];
            obj.process_noise = eye(numel(obj.state));
            obj.covariance = eye(numel(obj.state));
            obj.inov_covariance = zeros(numel(z),numel(z));
            obj.sizeZ = numel(z);
            obj.sizeX = numel(obj.state);
        end
        
        function z = makeMeasure(obj,x)
            %Predict measurement from state
            dof = numel(obj.model_handle.joints);
            obj.model_handle.position = x(1:dof);
            obj.model_handle.velocity = x(dof+1:dof*2);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
            
            %Generate measurements
            z = SensorMeasurement(obj.model_handle.sensors);
        end
        
        function x_new = stateUpdate(obj,x,u)
            %State Update
            
            dof = numel(obj.model_handle.joints);
                        
            x_new = x;
            %Update robot state by integrating
            x_new(1:dof) = x_new(1:dof)+...
                x(dof+1:dof*2)*u;
            
            %Forward Kinematics
            obj.model_handle.position = x_new(1:dof);
            obj.model_handle.velocity = x_new(dof+1:dof*2);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
        end
        
        function A = makeA(obj,x,u)
            dof = numel(obj.model_handle.joints);
            A = eye(numel(x),numel(x));
            A(1:dof,dof+1:dof*2) = A(1:dof,dof+1:dof*2)+eye(dof,dof)*u(1);
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