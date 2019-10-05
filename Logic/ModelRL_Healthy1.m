classdef ModelRL_Healthy1 < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, matPath)
            modelInstance = rlModelInstance_healthy1(0);
            modelInstance.loadModelFromModelSpecsNoSensor(xmlPath, matPath, 0);
            obj.model = modelInstance.model;
        end
    end
end