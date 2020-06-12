function ret = IOCRun(trialInfo, savePath, runParam)
    % Check function arguments and define default values if needed
    displayInfo = runParam.displayInfo;
    saveIntermediate = runParam.saveIntermediate;
    gamma = runParam.gamma;
    maxWinLen = runParam.maxWinLen;
    dimWeights = runParam.dimWeights;
    
    frameInds = 1:4000; %starting inds to run. leave empty to runa ll

    winSize = IOCInstanceNew.winSize;
    delta = IOCInstanceNew.delta;

    if ~exist('savePath', 'var')
        savePath = "";
    end
    
    % Return to calling function if no trial info is passed along
    if isempty(fieldnames(trialInfo))
        disp("No trial info was given");
        return
    end
    
    % Load to memory model and information needed for running IOC
    model = getModel(trialInfo.model, trialInfo.modelType, [], trialInfo.path);
%     model = ArmModelRL();
%     model.loadModelOnly();
% %     model.model.base = 'frame_6dof_root';
%     model.forwardKinematics();
    
    numWeights = length(trialInfo.candidateFeatures);
    numDofs = model.totalActuatedJoints;
    trueWeights = (trialInfo.weights)'/sum(trialInfo.weights);
    
    [q, dq, ddq, tau, states, control, trajT, trajU, trajX] = loadData(trialInfo, model);
    
    if 0 
        model.plotTrajectory(trajX);
    end

    % Initialization of threshold and condition variables
    dt = trajT(2) - trajT(1);

    % Create IOC instance
    ioc = IOCInstanceNew(model, dt);
    ioc.init(trialInfo);
    ioc.setFeatureNormalization(trajX, trajU);
    
    iocFeatures = ioc.calcFeatures(trajX, trajU);
    iocDynamics = ioc.calcDynamics(trajX, trajU);

    % Initialization of concatenated arrays summarizing error and gamma as functions of t
    if isempty(frameInds)
        frameInds = 1:size(trajT,2);
    end
    
    errorTraj = zeros(size(frameInds))';
    for indFrame = frameInds
        rankTraj = [];
        rankPassCodeTraj = [];
        
        H1 = [];
        H2 = [];
        
        for indSubFrame = indFrame:size(trajT,2)
			% start each individual window
			winInds = indFrame:indSubFrame;
            currLen = length(winInds);
            
            if currLen < winSize
                if displayInfo
                    fprintf("Obs frames %i:%i, width not sufficient. \n", ...
                        indFrame, indSubFrame);
                end
                continue;
            end
            
            % Read next observation
            x = trajX(winInds,:);
            u = trajU(winInds,:);

            % Generate new recovery matrix
            [H1, H2] = getRecoveryMatrix(ioc, H1, H2, x, u, dt);
            H = [H1 -H2];
            
            [checkHbool, checkHhave, checkHwant] = hPassFct1(H, dimWeights); 
            % checkHbool = 1;
            % checkHhave = 1;
            % checkHwant = 1;
            if ~checkHbool % run the inversion only when the window is properly sized
                if displayInfo
                    fprintf("Obs frames %i:%i, row to col ratio not sufficient. Have %u but want %u \n", ...
                        indFrame, indSubFrame, checkHhave, checkHwant);
                end
                continue;
            end
            
            % proceed with H matrix inversion
            Hhat = H/norm(H,'fro');
            
            % Compute weights using final recovery matrix
            [weights, ~] = computeWeights(Hhat, numWeights);
            % weightTraj(indFrame, :) = weights/sum(weights);
            
            % Compute error between true and recovered weights
            error = computeError(weights, trueWeights);
            errorTraj(indFrame,:) = error;
            
            % save rank of recovery matrix
            [rank, completed, rankPassCode] = validateCompletion(Hhat, gamma, delta, dimWeights);
            
            if currLen >= maxWinLen
                completed = 1;
            end
            
            rankTraj(currLen,:) = rank;
            rankPassCodeTraj(currLen, :) = rankPassCode;

            if displayInfo
                % Display information
                fprintf("Obs frames %i:%i of max %u/%u, rank: %0.2f/%0.2f, error: %0.3f, code: %u,%u,%u ", ...
                    indFrame, indSubFrame, currLen, maxWinLen, rank, gamma, error(1), ...
                    rankPassCode(1), rankPassCode(2), rankPassCode(3) );
                
                fprintf('Weights: ');
                for j = 1:numWeights
                    fprintf('%+0.3f,', weights(j));
                end
                fprintf('\n');
            end
            
            if completed
                % if it is done, save the final H matrix and weights, and clear out
                progressVar(indFrame).Hhat = Hhat;
                progressVar(indFrame).weights = weights'/sum(weights);
                progressVar(indFrame).rankTraj = rankTraj;
                progressVar(indFrame).winInds = winInds;
                break;
            end
        end

        % save either on the regular intervals, or at the last frame
        if (saveIntermediate > 0 && mod(indFrame, saveIntermediate) == 0) || indFrame == frameInds(end)
            outputVar.t = trajT;
            outputVar.q = q;
            outputVar.dq = dq;
            outputVar.tau = tau;
            outputVar.dt = dt;
            
            outputVar.lenDof = size(trajX, 2)/2;
            outputVar.lenState = size(trajX, 2);
            outputVar.lenControl = size(trajU, 2);
            
            outputVar.featureLabels = ioc.getFeatureListAsString();
            outputVar.features = iocFeatures;
            outputVar.dynamics = iocDynamics;
            
            outputVar.runParam = runParam;
            outputVar.trialInfo = trialInfo;
            outputVar.timeElapsed = toc;
            outputVar.timeElasedPerFrame = toc/length(indFrame);
            
            outputVar.frameInds = frameInds;
            outputVar.progress = progressVar;
            
            outputVar.minLenThres = checkHwant / outputVar.lenDof;
            outputVar.maxLenThres = maxWinLen;
            outputVar.minRankThres = gamma;
            
            % outputVar.errorTraj = errorTraj;
            % outputVar.rankTraj = rankTraj;
            % outputVar.weightTraj = weightTraj;
            % outputVar.completeTraj = completeTraj;
            % outputVar.rankPassCodeTraj = rankPassCodeTraj;
            
            % outputVar.segmentArray = segmentArray; % remove the first initialized value
            % outputVar.H1Log = H1Log;
            % outputVar.H2Log = H2Log;
            
            finalPath = sprintf("%s%s_%s_struct.mat", savePath, trialInfo.name, trialInfo.model);
            save(char(finalPath), 'outputVar');
        end
    end
end

function model = modifyModelWX(model)
    load('D:\aslab\projects\jf2lin\TROcopy\sub5.mat');
%     bodyPara = Para;

    bodyPara.lankle = 0.3700;
    bodyPara.lknee = 0.4000;
    bodyPara.lhip = 0.4500;
    bodyPara.mass = 65.7000;

    %legnth settings;
    l1=bodyPara.lankle;
    l2=bodyPara.lknee;
    l3=bodyPara.lhip;
    %mass settings
    m1=0.045*bodyPara.mass;      %0.09
    m2=0.146*bodyPara.mass;       %0.29
    m3=0.2985*bodyPara.mass;      %0.62
    %CoM position settings
    r1=(0.404)*l1;
    r2=(0.377)*l2;
    r3=(0.436)*l3;
    %moments of inertia settings
    I1=(0.28*l1)^2*m1;
    I2=(0.32*l2)^2*m2;             %0.35
    I3=(0.29*l3)^2*m3;             %0.30
    
    kinematicTransform(1).frameName = 'length_rknee_rankle';
    dynamicTransform(1).frameName = 'body_rknee_rankle';
    kinematicTransform(2).frameName = 'length_rhip_rknee';
    dynamicTransform(2).frameName = 'body_rhip_rknee';
    kinematicTransform(3).frameName = 'length_torso_rhip';
    dynamicTransform(3).frameName = 'body_torso_rhip';
    
%     kinematicTransform(4).frameName = 'length_rankle_rballfoot';
    
    kinematicTransform(1).t = eye(4);
    kinematicTransform(2).t = eye(4);
    kinematicTransform(3).t = eye(4);
%     kinematicTransform(4).t = eye(4);
    
    % determine which joint parameters use for defining upper-arm length
    % and second link inertia information
%     switch 'iit'
%         case 'iit'
            kinematicTransform(1).t(2, 4) = l1; % apply link length to the upper arm
            kinematicTransform(2).t(2, 4) = l2;
            kinematicTransform(3).t(2, 4) = l3;
%             
            dynamicTransform(1).m = m1;
            dynamicTransform(1).com = [0 r1 0]';
            dynamicTransform(1).I(1, 1) = I1;
            dynamicTransform(1).I(2, 2) = I1;
            dynamicTransform(1).I(3, 3) = 0;
            
            dynamicTransform(2).m = m2;
            dynamicTransform(2).com = [0 r2 0]';
            dynamicTransform(2).I(1, 1) = I2;
            dynamicTransform(2).I(2, 2) = I2;
            dynamicTransform(2).I(3, 3) = 0;
            
            dynamicTransform(3).m = m3;
            dynamicTransform(3).com = [0 r3 0]';
            dynamicTransform(3).I(1, 1) = I3;
            dynamicTransform(3).I(2, 2) = I3;
            dynamicTransform(3).I(3, 3) = 0;
            
%         case 'wj'
%             kinematicTransform(1).t(2, 4) = l1; % apply link length to the upper arm
%             kinematicTransform(2).t(2, 4) = l2;
%             kinematicTransform(3).t(2, 4) = l3;
%             
%             dynamicTransform(1).m = m1;
%             dynamicTransform(1).com = [0 -r1 0]';
%             dynamicTransform(1).I(1, 1) = I1;
%             dynamicTransform(1).I(2, 2) = 0;
%             dynamicTransform(1).I(3, 3) = I1;
%             
%             dynamicTransform(2).m = m2;
%             dynamicTransform(2).com = [0 -r2 0]';
%             dynamicTransform(2).I(1, 1) = I2;
%             dynamicTransform(2).I(2, 2) = 0;
%             dynamicTransform(2).I(3, 3) = I2;
%             
%             dynamicTransform(3).m = m3;
%             dynamicTransform(3).com = [0 -r3 0]';
%             dynamicTransform(3).I(1, 1) = I3;
%             dynamicTransform(3).I(2, 2) = 0;
%             dynamicTransform(3).I(3, 3) = I3;
%     end
    
    model.addKinDynTransformToModel(kinematicTransform, dynamicTransform);
    
    qOrig = q;
    dqOrig = dq;
    ddqOrig = ddq;
    
    q = [qOrig.ankle, qOrig.knee, qOrig.hip];
    dq = [dqOrig.ankle, dqOrig.knee, dqOrig.hip];
    ddq = [ddqOrig.ankle, ddqOrig.knee, ddqOrig.hip];
end

function [q, dq, ddq, tau, states, control, trajT, trajU, trajX] = loadData(trialInfo, model)
    % Load state, control and time trajectories to be analyzed.
    switch trialInfo.model
        case 'IITHalfBody'
            load(trialInfo.path);
            
            % keep only the joint angles corresponding
            qInds = [];
            allJointStr = {model.model.joints.name}';
             
            for indQ = 1:length(allJointStr)
                qInds(indQ) = find(ismember(saveVar.jointLabels, allJointStr{indQ}));
            end
           
            if ~isempty(trialInfo.runInds)
                tInds = trialInfo.runInds(1):trialInfo.runInds(2);
            else
                tInds = 1:length(saveVar.time);
            end

            time = saveVar.time(tInds);
            qRaw = saveVar.jointAngle.array(tInds, qInds);
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
%                 model.updateState(q(indTime, :), dq(indTime, :));
                tauRaw(indTime, :) = model.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
            end
            
%             tau = filter_dualpassBW(tauRaw, 0.04, 0, 5);
            tau = tauRaw;
            
%             states = [q dq];
            states = encodeState(q, dq);
            control = tau;
            
            trajT = time';
            trajU = control;
            trajX = states;
            
         otherwise
             inputPath = char(trialInfo.path);
             %inputPath = char(sprintf("%s%s", trialInfo.path, trialInfo.name));
             
            data = cell2mat(struct2cell(load(inputPath)));
            time = data(:,1);
            tau = data(:,2:numDofs+1);
            statesRaw = data(:, numDofs+2:end);
            
            q = statesRaw(:, 1:numDofs);
            dq = statesRaw(:, (numDofs+1):end);
            states = encodeState(q, dq);
            control = tau;
             
            % Approximate trajectories using splines
            trajT = linspace(0,1,1001);
            trajU = interp1(time, control, trajT,'spline');
            trajX = interp1(time, states, trajT,'spline');
    end    
end
