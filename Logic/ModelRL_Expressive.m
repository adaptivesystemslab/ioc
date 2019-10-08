classdef ModelRL_Expressive < ModelRL
    % ModelRL instance used to run IOC on single arm files of expressive
    % motion dataset
    
    properties
        modelType;
    end
    
    methods
        function obj = ModelRL_Expressive(type)        
            if nargin < 1
                obj.modelType = "ArmOnly";
            else
                obj.modelType = type;
            end
        end
        
        
        function loadModel(obj, xmlPath, matPath)
            if strcmp(obj.modelType, "ArmOnly")
                modelInstance = rlModelInstance_expressive_armOnly(0);
            else
                modelInstance = rlModelInstance_expressive_upperBody(0);
            end                      
            modelInstance.loadModelFromModelSpecsNoSensor(xmlPath, matPath, 1);
            obj.model = modelInstance.model;
        end
    end
end

