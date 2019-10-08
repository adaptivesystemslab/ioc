classdef ModelRL_Jumping_KW < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, trialInfo)
              obj.model = createJumpModel_ioc_2D_modForIOC(trialInfo.path, trialInfo.targNum, trialInfo.jumpNum, xmlPath);
%                     obj.modelJointNameRemap = {model_old.joints.name};
                    
%                     obj.model = model_old;
                    
%                     filepathModelInitPose = [obj.modelBaseFolder '/instance_jumping/model/JumpModel_Eul_inertia_2D_rightHalfBody.xml'];
%                     obj.model = createJumpModel_ioc_2D_modHalfBody(trialInfo.path, trialInfo.targNum, trialInfo.jumpNum, filepathModelInitPose);
%                     obj.model.forwardPosition();
%                     obj.model.base = 'rtoe0';
%                     obj.model.base = 'rframe1';
                    obj.model.forwardPosition();
        end
        
        function [q, dq, tau, trajT, trajU, trajX] = loadData(obj, trialInfo)
            load(trialInfo.path);
            
            targNum = trialInfo.targNum;
            jumpNum = trialInfo.jumpNum;
            
            param.jump.takeoffFrame = JA.TOFrame(jumpNum,targNum);
            param.jump.landFrame = JA.LandFrame(jumpNum,targNum);
            param.jump.locationLand = JA.locationLand(12*(targNum-1) + jumpNum);
            param.jump.grade = JA.jumpGrades(12*(targNum-1) + jumpNum);
            param.jump.modelLinks = JA.modelLinks;
            param.jump.world2base = squeeze(JA.world2base((12*(targNum-1) + jumpNum),:,:));
            param.jump.bad2d = JA.bad2D(jumpNum,targNum);
            
                    % Crop out initial and final calibration motions
%                     takeoffFrames = 200; % 1 second before takeoff frame ...
%                     landFrames = 300; % ... to 1.5 seconds after takeoff frame (~ 1 second after landing)
%                     framesToUse = (param.jump.takeoffFrame-takeoffFrames):(param.jump.takeoffFrame+landFrames);
                    
            takeoffFrames = 0; % 1 second before takeoff frame ...
            landFrames = 0; % ... to 1.5 seconds after takeoff frame (~ 1 second after landing)
            framesToUse = (param.jump.takeoffFrame-takeoffFrames):(param.jump.landFrame+landFrames);

            if(numel(framesToUse) < size(JA.targ(targNum).jump(jumpNum).data,1))
                fullDataAngles = JA.targ(targNum).jump(jumpNum).data(framesToUse,:);
            else % jump recording stops sooner than "landFrames" after TOFrame
                fullDataAngles = JA.targ(targNum).jump(jumpNum).data( framesToUse(1):end ,:);
                fullDataAngles = [fullDataAngles; repmat(fullDataAngles(end,:),(numel(framesToUse) - size(fullDataAngles,1)),1)]; % repeat last joint angle measurement for remainder of frames
            end

            % keep only a subset of the joint angles
%             qInds = [];
%             allJointStr = {obj.model.joints.name}';
%             for indQ = 1:length(allJointStr)
%                 qInds(indQ) = find(ismember(obj.modelJointNameRemap, allJointStr{indQ}));
%             end
            qInds = 1:length(obj.model.joints);
            
            % also, negate the following joints since they're past
            % the flip
            qFlip = fullDataAngles;
            %                     jointsToFlip = {'rankle_jDorsiflexion', 'rknee_jExtension', 'rhip_jFlexion'};
            % %                     jointsToFlip = {'rankle_jDorsiflexion', 'rknee_jExtension', 'rhip_jFlexion', 'back_jFB', 'rjoint1'};
            %                     for indQ = 1:length(jointsToFlip)
            %                         qIndsFlip(indQ) = find(ismember(model.modelJointNameRemap, jointsToFlip{indQ}));
            %                         qFlip(:, qIndsFlip(indQ)) = -fullDataAngles(:, qIndsFlip(indQ));
            %                     end

            dt = 0.005;
            time = dt*(0:(size(qFlip, 1)-1));
            qRaw = qFlip(:, qInds);
            q = filter_dualpassBW(qRaw, 0.04, 0, 5);

            dqRaw = calcDerivVert(q, dt);
            dq = filter_dualpassBW(dqRaw, 0.04, 0, 5);
            %             dq = dqRaw;

            % don't filter ddq and tau to keep
            ddqRaw = calcDerivVert(dq, dt);
            %             ddq = filter_dualpassBW(ddqRaw, 0.04, 0, 5);
            ddq = ddqRaw;

            tauRaw = zeros(size(q));
            for indTime = 1:length(time) % recalc torque given redistributed masses
                tauRaw(indTime, :) = obj.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
            end

            %             tau = filter_dualpassBW(tauRaw, 0.04, 0, 5);
            tau = tauRaw;

            %             states = [q dq];
            states = encodeState(q, dq);
            control = tau;

            trajT = time';
            trajU = control;
            trajX = states;
        end
    end
end