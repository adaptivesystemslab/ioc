classdef EKF_Q_DQ_DDQ < EKF_Resize
    
    properties
        %Robot Model Handle
        %model_handle = [];
    end
    
    methods
        function obj = EKF_Q_DQ_DDQ(model)
            
            %Call super
            obj = obj@EKF_Resize(model);
            %obj.model_handle = model;
            %Our observations are the sensors attached to the model
            z = vertcat(model.sensors.measurement);
            obj.observation_noise = eye(numel(z));
            %Our state is model q dq ddq
            obj.state = [...
                model.position ; ...
                model.velocity ; ...
                model.acceleration ; ...
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
            obj.model_handle.acceleration = x(dof*2+1:dof*3);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
            obj.model_handle.forwardAcceleration();
            
            %Generate measurements
            z = SensorMeasurement(obj.model_handle.sensors);
        end
        
        function x_new = stateUpdate(obj,x,u)
            %State Update
            
            dof = numel(obj.model_handle.joints);
                        
            x_new = x;
            %Update robot state by integrating
            x_new(1:dof) = x_new(1:dof)+...
                x(dof+1:dof*2)*u + x(dof*2+1:dof*3)*u^2/2;
            x_new(dof+1:dof*2) = x_new(dof+1:dof*2) + x(dof*2+1:dof*3)*u;
            
            %Forward Kinematics
            obj.model_handle.position = x_new(1:dof);
            obj.model_handle.velocity = x_new(dof+1:dof*2);
            obj.model_handle.acceleration = x_new(dof*2+1:dof*3);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
            obj.model_handle.forwardAcceleration();
        end
        
        function A = makeA(obj,x,u)
            dof = numel(obj.model_handle.joints);
            A = eye(numel(x),numel(x));
            A(1:dof*2,dof+1:dof*3) = A(1:dof*2,dof+1:dof*3)+eye(dof*2,dof*2)*u(1);
            A(1:dof,dof*2+1:dof*3) = A(1:dof,dof*2+1:dof*3)+eye(dof,dof)*u(1)*u(1)/2;
        end
         
        function H = makeH(obj,x)
            %Generate Sensor Jacobian Matrix based on sensors attached to
            %the model
            obj.model_handle.calculateSensorJacobians();
            H = vertcat(obj.model_handle.sensors.obsJacobian);
        end
    end
end