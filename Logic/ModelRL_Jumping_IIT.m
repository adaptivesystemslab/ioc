classdef ModelRL_Jumping_IIT < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, matPath)
            modelInstance = rlModelInstance_jumping(0);
            modelInstance.loadModelFromModelSpecsNoSensor(xmlPath, matPath, 1);
            obj.model = modelInstance.model;
        end
    end
end