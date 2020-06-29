classdef EKF_Constraint < handle & matlab.mixin.Heterogeneous
    %Base class for a constraint definition used by Jacobian constrained
    %EKF. Each constraint is really a pair of sensors attached to the
    %models without decorators.
    
    properties
        models = rlCModel.empty();
        sensors = SensorCore.empty();
        A = [];
        b = [];
    end
    
    
    methods
        function obj = EKF_Constraint(name,model1,frame1,T1,model2,frame2,T2)
            %Constructor creates a pair of sensors and attaches them to the models
            
            %Save model references
            obj.models = [model1;model2];
            %Create the decoratorless sensors
            obj.sensors(1) = SensorCore([name '_A']);
            obj.sensors(2) = SensorCore([name '_B']);
            %Attach the sensors to the models
            model1.addSensor(obj.sensors(1),frame1,T1);
            model2.addSensor(obj.sensors(2),frame2,T2);
            %Run forward position 
            model1.forwardPosition();
            model2.forwardPosition();
        end
        
        function value = get.A(obj)
           value = obj.get_A_value(); 
        end
        
        function value = get.b(obj)
           value = obj.get_b_value(); 
        end
        
    end
    
    methods (Access = protected)
        function value = get_A_value(obj)
            value = [];
        end
        
        function value = get_b_value(obj)
           value = []; 
        end
    end
end