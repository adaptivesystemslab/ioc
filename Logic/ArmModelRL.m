classdef ArmModelRL < handle

    properties
        model;
        model_old;
        modelJointNameRemap;
%         modelBaseFolder = '../../kalmanfilter/ik_framework';
        modelBaseFolder = '../Libraries/rl/ik_framework';
        
        totalActuatedJoints;
    end
    
    methods
        function obj = ArmModelRL()        
            
        end
        
        function initModel(obj)
            obj.totalActuatedJoints = length(obj.model.joints);
        end
        
        function createModel(obj, jointParams)
            % create model according to Peter Corke Toolbox
            joints = ModelJoint(jointParams(1));
            for i = 2:length(jointParams)
                joints(i) = ModelJoint(jointParams(i));
            end
            
            switch length(joints)
                case 3
                    modelStr = '../RLModels/ioc_v4_3dof.xml';
                    obj.model = rlCModel(modelStr);
                    
                case 5
                    modelStr = '../RLModels/ioc_v4_4dof.xml';
                    obj.model = rlCModel(modelStr);
            end
            
            switch length(joints)
                case 3
                    [kinematicTransform, dynamicTransform] = obj.setup3DofModel(joints);
                    obj.addKinDynTransformToModel(kinematicTransform, dynamicTransform);
                    
                case 5
                    [kinematicTransform, dynamicTransform] = obj.setup4DofModel(joints);
                    obj.addKinDynTransformToModel(kinematicTransform, dynamicTransform);
            end
        end
        
        function [kinematicTransform, dynamicTransform] = setup3DofModel(obj, joints)
              % apply link offsets to the model
            kinematicTransform(1).frameName = 'length_rshoulder_relbow';
            dynamicTransform(1).frameName = 'body_rshoulder_relbow';
            
            kinematicTransform(2).frameName = 'length_relbow_rwrist';
            dynamicTransform(2).frameName = 'body_relbow_rwrist';
            
            kinematicTransform(1).t = eye(4);
            kinematicTransform(2).t = eye(4);
            
            % determine which joint parameters use for defining upper-arm length
            % and second link inertia information
            kinematicTransform(1).t(1, 4) = joints(3).length; % apply link length to the upper arm
            %kinematicTransform(2).t(1, 4) = obj.joints(3).length; % apply link length to the fore arm
            
            dynamicTransform(1).m = joints(2).mass;
            dynamicTransform(1).com = joints(2).com;
            dynamicTransform(1).I = joints(2).inertiaTensor;
            
            dynamicTransform(2).m = joints(3).mass;
            dynamicTransform(2).com = joints(3).com;
            dynamicTransform(2).I = joints(3).inertiaTensor;
        end
        
        function [kinematicTransform, dynamicTransform] = setup4DofModel(obj, joints)
              % apply link offsets to the model
            kinematicTransform(1).frameName = 'length_rshoulder_relbow';
            dynamicTransform(1).frameName = 'body_rshoulder_relbow';
            
            kinematicTransform(2).frameName = 'length_relbow_rwrist';
            dynamicTransform(2).frameName = 'body_relbow_rwrist';
            
            kinematicTransform(1).t = eye(4);
            kinematicTransform(2).t = eye(4);
            
            % determine which joint parameters use for defining upper-arm length
            % and second link inertia information
            kinematicTransform(1).t(3, 4) = joints(3).offset; % apply link length to the upper arm
            kinematicTransform(2).t(3, 4) = joints(5).offset;
            
            dynamicTransform(1).m = joints(3).mass;
            dynamicTransform(1).com = (-1*joints(3).com)';
            dynamicTransform(1).I = joints(3).inertiaTensor;
            
            dynamicTransform(2).m = joints(5).mass;
            dynamicTransform(2).com = (-1*joints(5).com)';
            dynamicTransform(2).I = joints(5).inertiaTensor;
        end
        
        function loadJumpingModel(obj, filepathSourceMat, trialInfo)
            switch trialInfo.model
                case 'Jumping'
                    % load original IIT-based RL form
                    filepathModelInitPose = [obj.modelBaseFolder '/instance_jumping/model/JumpModel_IIT.xml'];
                    modelInstance = rlModelInstance_jumping(0);
                    modelInstance.loadModelFromModelSpecsNoSensor(filepathModelInitPose, filepathSourceMat, 1);
                    obj.model = modelInstance.model;
                    
                case 'Jumping2D'
                    % load Kevin's modified form
                    filepathModelInitPose = [obj.modelBaseFolder '/instance_jumping/model/JumpModel_Eul_inertia_2D.xml'];
                    obj.model_old = createJumpModel_ioc_2D_modForIOC(trialInfo.path, trialInfo.targNum, trialInfo.jumpNum, filepathModelInitPose);
                    obj.modelJointNameRemap = {obj.model_old.joints.name};
                    
                    filepathModelInitPose = [obj.modelBaseFolder '/instance_jumping/model/JumpModel_Eul_inertia_2D_rightHalfBody.xml'];
                    obj.model = createJumpModel_ioc_2D_modHalfBody(trialInfo.path, trialInfo.targNum, trialInfo.jumpNum, filepathModelInitPose);
                    obj.model.forwardPosition();
                    obj.model.base = 'rtoe0';
%                     obj.model.base = 'rframe1';
                    obj.model.forwardPosition();
            end
        end
        
        function loadIITModel(obj, filepathSourceMat, trialInfo)
            switch trialInfo.model
                case 'IIT_17DOF'
                    filepathModelInitPose = [obj.modelBaseFolder '/instance_iit/model/iit_v10_fixedbase_right.xml']; % 3D, 17dof half body model
            
                case 'IIT_7DOF'
                    filepathModelInitPose = [obj.modelBaseFolder '/instance_iit/model/iit_v10_fixedbase_right_sag.xml']; % sag 2D, 7dof half body model
             
                case 'IIT_3DOF'
                    filepathModelInitPose = [obj.modelBaseFolder '/instance_iit/model/iit_v10_fixedbase_rightleg_sag.xml']; % sag 2D, 3dof half body model
            end
            
%             filepathModelInitPose = '../../kalmanfilter/ik_framework/instance_iit/model/iit_v10_fixedbase_rightlegonly.xml'; % sag 2D, 3dof leg only
     
            modelInstance = rlModelInstance_iit(0);
            modelInstance.loadModelFromModelSpecsNoSensor(filepathModelInitPose, filepathSourceMat);
            obj.model = modelInstance.model;
        end
        
        function loadModelOnly(obj)
            filepathModelInitPose = '../../kalmanfilter/ik_framework/instance_iit/model/iit_v10_fixedbaseRightAnkle.xml';
            obj.model = rlCModel(filepathModelInitPose);
            
            if 0
                mdl = obj.model;
                vis = rlVisualizer('vis',640,480);
                mdl.forwardPosition();
                vis.addModel(mdl);
                vis.update();
                
                mdl.position = [0 0 0];
                mdl.forwardPosition();
                vis.addMarker('x-axis', [1 0 0], [0.2 0.2 0.2 1]);
                vis.addMarker('y-axis', [0 1 0], [0.2 0.2 0.2 1]);
                vis.addMarker('z-axis', [0 0 1], [0.2 0.2 0.2 1]);
                for i = 1:length(mdl.joints)
                    vis.addMarker(['j_' mdl.joints(i).name], ...
                        mdl.joints(i).t(1:3, 4), [0.4 0.4 0.8 1]);
                end
                for i = 1:length(mdl.bodies)
                    vis.addMarker(['EF_' mdl.bodies(i).name], ...
                        mdl.bodies(i).t(1:3, 4), [0.4 0.8 0 1]);
                    vis.addMarker(['COM_' mdl.bodies(i).name], ...
                        mdl.bodies(i).t(1:3, 4) + mdl.bodies(i).com, [0.8 0 0.4 1]);
                end
                vis.update();
            end
        end
        
        function loadAndSetupJumping(obj, filepathSourceMat, trialInfo)
            % load the model
            obj.loadJumpingModel(filepathSourceMat, trialInfo);
%             obj.model.base = 'frame_6dof_root';
            obj.forwardKinematics();

%            JumpModel_IIT.xml 
        end
        
        function loadAndSetupIIT(obj, filepathSourceMat, trialInfo)
            % present patient variables
            height = 1.70;
            weight = 62;
            gender = 'f';
            
            % load the model
            obj.loadIITModel(filepathSourceMat, trialInfo);
            obj.model.base = 'frame_6dof_root';
            obj.forwardKinematics();
            
            allBodyNames = {obj.model.bodies.name};

            % add head mass to torso
            fprintf('Merging headneck mass into torso\n');
            dumasFrameStr = 'head&neck';
            mass =            lookupTableDumas('mass',     dumasFrameStr, gender, [], [])*weight;
            bodyInd = find(ismember(allBodyNames, 'body_l5s1_t1c7'));
            obj.model.bodies(bodyInd).m = obj.model.bodies(bodyInd).m + mass;
            
            % double right arm and leg to compensate for missing left arm
            % and leg
            fprintf('Doubling right leg and arm segment masses to compensate for truncated model\n');
            bodiesToDouble = {'body_rhip_rknee', 'body_rknee_rankle', 'body_rankle_rballfoot', ...
                'body_c7rshoulder_rshoulder', 'body_rshoulder_relbow', 'body_relbow_rwrist', 'body_rwrist_rhand'};
            for i = 1:length(bodiesToDouble)
                bodyInd = find(ismember(allBodyNames, bodiesToDouble{i}));
                obj.model.bodies(bodyInd).m = obj.model.bodies(bodyInd).m*2;
            end
        end
        
        function addEndEffectors(obj, frameStr)
            % add a sensor to the end effector to denote where it is
            for i = 1:length(frameStr)
                M1 = SensorCore(['endEff' num2str(i) '_' frameStr{i}]);
                obj.model.addSensor(M1, frameStr{i}, eye(4));
            end
            
            clear M1;
        end
        
        function qVec = getState(obj)
            q = obj.model.position;
            dq = obj.model.velocity;
            qVec = [q dq]';
        end
        
        function updateState(obj, newQ, newDQ)
            obj.model.position = newQ;
            obj.model.velocity = newDQ;
            
            obj.forwardKinematics();
        end
        
        function updateAccelerations(obj, ddq)
            obj.model.acceleration = ddq;
        end
        
        function ddq = getAcceleration(obj)
            ddq = obj.model.acceleration;
        end
        
        function forwardKinematics(obj)
            obj.model.forwardPosition();
            obj.model.forwardVelocity();
            obj.model.forwardAcceleration();
        end
        
        function x = getEndEffectorPosition(obj, i)
            x = obj.model.sensors(i).transform(1:3, 4);
        end
        
        function dx = getEndEffectorVelocity(obj, i)
            obj.forwardKinematics();
            obj.model.calculateSensorJacobians();
            J = obj.model.sensors(i).baseJacobian;
            
            stateVec = obj.getState();
            dq = stateVec(2, :);
            dxTemp = J*dq';
            dx = dxTemp(4:6)';
        end
        
        function ddx = getEndEffectorAcceleration(obj, ddq, i)
            obj.forwardKinematics();
            obj.model.calculateSensorJacobians();
            J = obj.model.sensors(i).baseJacobian;
            dJ = obj.model.sensors(i).baseJacobianDerivative;
            
            stateVec = obj.getState();
            dq = stateVec(2, :);
            ddxTemp = dJ*dq' + J*ddq';
            ddx = ddxTemp(4:6)';
        end
        
        function torques = inverseDynamics(obj, accelerations)
            stateVec = obj.getState();
            obj.model.position = stateVec(1,:)';
            obj.model.velocity = stateVec(2,:)';
            obj.model.acceleration = accelerations;
            
            obj.forwardKinematics();
            
            obj.model.inverseDynamics();
            torques = obj.model.torque';
        end
        
        function torques = inverseDynamicsQDqDdq(obj, q, dq, ddq)
            obj.model.position = q;
            obj.model.velocity = dq;
            obj.model.acceleration = ddq;
            
            obj.forwardKinematics();
            obj.model.inverseDynamics();
            
            torques = obj.model.torque';
        end
        
        function ddq = forwardDynamicsQDqTau(obj, q, dq, tau)
            obj.model.position = q;
            obj.model.velocity = dq;
            obj.model.torque = tau;
            
            obj.forwardKinematics();
            obj.model.forwardDynamics();
            
            ddq = obj.model.acceleration';
        end
        
        function ddqVec = forwardDynamics(obj, torques, doUpdate)
            if ~exist('doUpdate', 'var')
                doUpdate = 0;
            end
            
            % save the old state in case doUpdate is not checked
            if ~doUpdate
                oldPos = obj.model.position;
                oldVel = obj.model.velocity;
                oldAcc = obj.model.acceleration;
            end
            
            obj.model.torque = torques;
            obj.model.forwardDynamics();
            
            ddqVec = obj.model.acceleration';
            
            if ~doUpdate
                obj.model.position = oldPos;
                obj.model.velocity = oldVel;
                obj.model.acceleration = oldAcc;
                
                obj.forwardKinematics();
            end
        end
        
        function inertialMatrix = getInertialMatrix(obj)
            stateVec = obj.getState();
            obj.model.position = stateVec(1,:);
            obj.model.velocity = stateVec(2,:);
            
            obj.forwardKinematics();
            
            obj.model.calculateMassMatrix();
            
            inertialMatrix = obj.model.M;
        end
        
        function addKinDynTransformToModel(obj, kinematicTransform, dynamicTransform)
            for i = 1:length(kinematicTransform)
                if ~isempty(kinematicTransform(i).t)
                    obj.setKinematicTransform(kinematicTransform(i));
                end
            end
            
            for i = 1:length(dynamicTransform)
                if ~isempty(dynamicTransform(i).m)
                    obj.setDynamicTransform(dynamicTransform(i));
                end
            end
        end
        
        function setKinematicTransform(obj, kinematicTransform)
            allFrameNames = {obj.model.transforms.name};
            indFrame = find(ismember(allFrameNames, kinematicTransform.frameName) == 1);
            
            if ~isempty(indFrame)
                obj.model.transforms(indFrame).t = kinematicTransform.t;
            else
                fprintf('Frame not found: %s\n', kinematicTransform.frameName);
             end
        end
        
        function setDynamicTransform(obj, dynamicTransform)
            allFrameNames = {obj.model.bodies.name};
            indFrame = find(ismember(allFrameNames, dynamicTransform.frameName) == 1);
            
            if ~isempty(indFrame)
                obj.model.bodies(indFrame).m   = dynamicTransform.m;
                obj.model.bodies(indFrame).com = dynamicTransform.com;
                obj.model.bodies(indFrame).I   = dynamicTransform.I;
            else
                fprintf('Frame not found: %s\n', dynamicTransform.frameName);
            end
        end
        
        function result = getDerivativesRelativeToState(obj, control)
            % Computes finite difference derivatives of the forward
            % dynamic equations wrt each state variable
            
            % Dertivatives are computed numerically using central finite
            % differences
                        
            delta = 5e-2;
            % dofs = length(obj.joints);
            dofs = obj.totalActuatedJoints;
            result = zeros(dofs,dofs*2);
            
            currentState = obj.getState();
                        
            for i = 0:dofs-1
                changes = zeros(1,dofs);
                changes(i+1) = delta;
                                       
                % Define right and left states for joint i
                leftStateQ = currentState(1,:) - changes;
                rightStateQ = currentState(1,:) + changes;
                leftStateDQ = currentState(2,:) - changes;
                rightStateDQ = currentState(2,:) + changes;

                % Calculate finite difference wrt angle state variable
                obj.updateState(rightStateQ, currentState(2,:));
                ddq = obj.forwardDynamics(control);
                result(:, 2*i+1) = ddq';
                    
                obj.updateState(leftStateQ, currentState(2,:));
                ddq = obj.forwardDynamics(control);
                result(:, 2*i+1) = (result(:, 2*i+1) - ddq')/(2*delta);
                    
                 % Calculate finite difference wrt velocity state
                 % variable
                 obj.updateState(currentState(1,:), rightStateDQ);
                 ddq = obj.forwardDynamics(control);
                 result(:, 2*i+2) = ddq';
                    
                 obj.updateState(currentState(1,:), leftStateDQ);
                 ddq = obj.forwardDynamics(control);
                 result(:, 2*i+2) = (result(:, 2*i+2) - ddq')/(2*delta);
            end
            
            obj.updateState(currentState(1,:), currentState(2,:));
        end
        
        function plotTrajectory(obj, q)
            vis = rlVisualizer('vis',640,480);
            obj.model.forwardPosition();
            vis.addModel(obj.model);
            vis.update();
            
            for i = 1:size(q, 1)
                obj.model.position = q(i, :);
                obj.model.forwardPosition();
                
                pause(0.001);
                vis.update();
            end
        end
    end   
 end