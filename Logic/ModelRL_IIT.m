classdef ModelRL_IIT < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, matPath)
            modelInstance = rlModelInstance_iit(0);
            modelInstance.loadModelFromModelSpecsNoSensor(xmlPath, matPath);
            obj.model = modelInstance.model;
            
            % switching the base to the ankle
            obj.model.base = 'frame_6dof_root';
            obj.forwardKinematics();

            % add head mass to torso since this model only loads half the body
            allBodyNames = {obj.model.bodies.name};
            fprintf('Merging headneck mass into torso\n');
            dumasFrameStr = 'head&neck';
            mass =            lookupTableDumas('mass',     dumasFrameStr, modelInstance.gender, [], [])*modelInstance.body_weight;
            bodyInd = find(ismember(allBodyNames, 'body_l5s1_t1c7'));
            obj.model.bodies(bodyInd).m = obj.model.bodies(bodyInd).m + mass;
            
            % double right arm and leg to compensate for missing left arm and leg
            fprintf('Doubling right leg and arm segment masses to compensate for truncated model\n');
            bodiesToDouble = {'body_rhip_rknee', 'body_rknee_rankle', 'body_rankle_rballfoot', ...
                'body_c7rshoulder_rshoulder', 'body_rshoulder_relbow', 'body_relbow_rwrist', 'body_rwrist_rhand'};
            for i = 1:length(bodiesToDouble)
                bodyInd = find(ismember(allBodyNames, bodiesToDouble{i}));
                obj.model.bodies(bodyInd).m = obj.model.bodies(bodyInd).m*2;
            end
        end
    end
end