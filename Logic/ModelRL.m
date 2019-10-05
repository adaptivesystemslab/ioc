classdef ModelRL < handle
    properties
        model;

%         modelJointNameRemap;
%         modelBaseFolder = '../Libraries/rl/ik_framework';
        
%         totalActuatedJoints;
    end
    
    methods
        function obj = ModelRL()        
            
        end
        
        function numDof = getModelDof(obj)
            numDof = length(obj.model.joints);
        end
        
        function [q, dq, ddq, tau, trajT, trajU, trajX, traj] = loadData(obj, trialInfo)
            load(trialInfo.path);
            
            % keep only the joint angles corresponding
            qInds = [];
            allJointStr = {obj.model.joints.name}';
            
            for indQ = 1:length(allJointStr)
                qInds(indQ) = find(ismember(saveVar.jointLabels, allJointStr{indQ}));
            end
            
            time = saveVar.time;
            qRaw = saveVar.jointAngle.array(:, qInds);
            q = filter_dualpassBW(qRaw, 0.04, 0, 5);
            
            dqRaw = calcDerivVert(q, saveVar.dt);
            dq = filter_dualpassBW(dqRaw, 0.04, 0, 5);
%             dq = dqRaw;
            
            % don't filter ddq and tau to keep
            ddqRaw = calcDerivVert(dq, saveVar.dt);
%             ddq = filter_dualpassBW(ddqRaw, 0.04, 0, 5);
            ddq = ddqRaw;
            
            tauRaw = zeros(size(q));
            for indTime = 1:length(time) % recalc torque given redistributed masses
                tauRaw(indTime, :) = obj.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
            end
            
%             tau = filter_dualpassBW(tauRaw, 0.04, 0, 5);
            tau = tauRaw;
            
            states = encodeState(q, dq);
            control = tau;
            
            trajT = time';
            trajU = control;
            trajX = states;
            
            traj.q=q;
            traj.dq=dq;
            traj.ddq=ddq;
            traj.tau=tau;
            traj.states=states;
            traj.control=control;
            traj.time=time';
            traj.trajT=trajT;
            traj.trajU=trajU;
            traj.trajX=trajX;
            traj.frameInds=1:length(time);            
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
            
            obj.forwardKinematics();
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
            
            dq = obj.model.velocity;
            dxTemp = J*dq';
            dx = dxTemp(4:6)';
        end
        
        function ddx = getEndEffectorAcceleration(obj, ddq, i)
            obj.forwardKinematics();
            obj.model.calculateSensorJacobians();
            J = obj.model.sensors(i).baseJacobian;
            dJ = obj.model.sensors(i).baseJacobianDerivative;
            
            dq = obj.model.velocity;
            ddxTemp = dJ*dq + J*ddq';
            ddx = ddxTemp(4:6)';
        end
        
%         function torques = inverseDynamics(obj, accelerations)
%             obj.model.acceleration = accelerations;
%             obj.forwardKinematics();
%             
%             obj.model.inverseDynamics();
%             torques = obj.model.torque';
%         end
        
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
        
%         function ddqVec = forwardDynamics(obj, torques, doUpdate)
%             if ~exist('doUpdate', 'var')
%                 doUpdate = 0;
%             end
%             
%             % save the old state in case doUpdate is not checked
%             if ~doUpdate
%                 oldPos = obj.model.position;
%                 oldVel = obj.model.velocity;
%                 oldAcc = obj.model.acceleration;
%             end
%             
%             obj.model.torque = torques;
%             obj.model.forwardDynamics();
%             
%             ddqVec = obj.model.acceleration';
%             
%             if ~doUpdate
%                 obj.model.position = oldPos;
%                 obj.model.velocity = oldVel;
%                 obj.model.acceleration = oldAcc;
%                 
%                 obj.forwardKinematics();
%             end
%         end
        
        function inertialMatrix = getInertialMatrix(obj)
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
        
%         function result = getDerivativesRelativeToState(obj, control)
%             % Computes finite difference derivatives of the forward
%             % dynamic equations wrt each state variable
%             
%             % Dertivatives are computed numerically using central finite
%             % differences
%                         
%             delta = 5e-2;
%             % dofs = length(obj.joints);
%             dofs = obj.totalActuatedJoints;
%             result = zeros(dofs,dofs*2);
%             
%             currentState = obj.getState();
%                         
%             for i = 0:dofs-1
%                 changes = zeros(1,dofs);
%                 changes(i+1) = delta;
%                                        
%                 % Define right and left states for joint i
%                 leftStateQ = currentState(1,:) - changes;
%                 rightStateQ = currentState(1,:) + changes;
%                 leftStateDQ = currentState(2,:) - changes;
%                 rightStateDQ = currentState(2,:) + changes;
% 
%                 % Calculate finite difference wrt angle state variable
%                 obj.updateState(rightStateQ, currentState(2,:));
%                 ddq = obj.forwardDynamics(control);
%                 result(:, 2*i+1) = ddq';
%                     
%                 obj.updateState(leftStateQ, currentState(2,:));
%                 ddq = obj.forwardDynamics(control);
%                 result(:, 2*i+1) = (result(:, 2*i+1) - ddq')/(2*delta);
%                     
%                  % Calculate finite difference wrt velocity state
%                  % variable
%                  obj.updateState(currentState(1,:), rightStateDQ);
%                  ddq = obj.forwardDynamics(control);
%                  result(:, 2*i+2) = ddq';
%                     
%                  obj.updateState(currentState(1,:), leftStateDQ);
%                  ddq = obj.forwardDynamics(control);
%                  result(:, 2*i+2) = (result(:, 2*i+2) - ddq')/(2*delta);
%             end
%             
%             obj.updateState(currentState(1,:), currentState(2,:));
%         end
        
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