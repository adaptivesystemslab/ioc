classdef ModelRL_IIT < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, matPath)
            modelInstance = rlModelInstance_iit(0);
            modelInstance.loadModelFromModelSpecsNoSensor(xmlPath, matPath);
            obj.model = modelInstance.model;
        end
    end
end