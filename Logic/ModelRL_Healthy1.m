classdef ModelRL_Healthy1 < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, matPath)
            modelInstance = rlModelInstance_healthy1(0);
            saveVar = modelInstance.loadModelFromModelSpecsNoSensor(xmlPath, matPath, 0);
            obj.model = modelInstance.model;
            
            obj.model.base = saveVar.baseFrame;
            obj.model.forwardPosition();
        end
    end
end