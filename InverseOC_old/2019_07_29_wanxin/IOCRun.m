function ret = IOCRun(trialInfo, savePath, runParam)
    % Check function arguments and define default values if needed
    % displayInfo = runParam.displayInfo;
    % saveIntermediate = runParam.saveIntermediate;
    % gamma = runParam.gamma;
    % maxWinLen = runParam.maxWinLen;
    % dimWeights = runParam.dimWeights;
    
    frameInds = 1:1000; %starting inds to run. leave empty to run all
    
    runParam.hWinAdvFlag = 'forward'; % forward, backward, centre

    runParam.winSize = IOCInstanceNew.winSize;
    runParam.delta = IOCInstanceNew.delta;

    if ~exist('savePath', 'var')
        savePath = "";
    end
    
    % Return to calling function if no trial info is passed along
    if isempty(fieldnames(trialInfo))
        disp("No trial info was given");
        return
    end
    
    % Load to memory model and information needed for running IOC
    if 1 % wanxin modification trial
        model = ArmModelRL();
        model.loadModelOnly();
        [q, dq, ddq, tau, states, control, trajT, trajU, trajX, model] = modifyModelWX(model);
    else
        model = getModel(trialInfo.model, trialInfo.modelType, [], trialInfo.path);
        [q, dq, ddq, tau, states, control, trajT, trajU, trajX] = loadData(trialInfo, model);
    end
    
    runParam.numWeights = length(trialInfo.candidateFeatures);
    runParam.numDofs = model.totalActuatedJoints;
    runParam.trueWeights = (trialInfo.weights)'/sum(trialInfo.weights);

    % Initialization of threshold and condition variables
    dt = trajT(2) - trajT(1);

    % Create IOC instance
    ioc = IOCInstanceNew(model, dt);
    ioc.init(trialInfo);
%     ioc.setFeatureNormalization(trajX, trajU);
    
    iocFeatures = ioc.calcFeatures(trajX, trajU);
    iocDynamics = ioc.calcDynamics(trajX, trajU);

    % Initialization of concatenated arrays summarizing error and gamma as functions of t
    if isempty(frameInds)
        frameInds = 1:size(trajT,2);
    end
    
    runParam.lenDof = size(trajX, 2)/2;
    runParam.hWant = (size(trajX, 2) + runParam.numWeights) * runParam.dimWeights;

    % precalc singular H1 and H2 matrix so 
    trajH1 = [];
    trajH2 = [];
    for indFrame = frameInds
        switch runParam.displayInfo
            case 'verbose'
                fprintf("Obs frames %i, precalculating H1/H2... \n", ...
                    indFrame);
        end
        
        [H, H1, H2] = assembleHMatrix(trajT, trajX, trajU, dt, ioc, [], [], indFrame, runParam);
        
        trajH1(:, :, indFrame) = H1;
        trajH2(:, :, indFrame) = H2;
    end
    
    % assemble and assess H
    for indFrame = frameInds
        progressVar(indFrame) = calcWinLenAndH(trajT, trajX, trajU, dt, ioc, indFrame - 1, trajH1, trajH2, runParam);
        
        % save either on the regular intervals, or at the last frame
        if (runParam.saveIntermediate > 0 && mod(indFrame, runParam.saveIntermediate) == 0) || indFrame == frameInds(end)
            outputVar.t = trajT;
            outputVar.q = q;
            outputVar.dq = dq;
            outputVar.tau = tau;
            outputVar.dt = dt;
            
            outputVar.lenDof = runParam.lenDof;
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
            
            outputVar.minLenThres = runParam.hWant / outputVar.lenDof;
            outputVar.maxLenThres = runParam.maxWinLen;
            outputVar.minRankThres = runParam.gamma;
            
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

function [q, dq, ddq, tau, states, control, trajT, trajU, trajX, model] = modifyModelWX(model)
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
        
    
    model.addKinDynTransformToModel(kinematicTransform, dynamicTransform);
    
    qOrig = q;
    dqOrig = dq;
    ddqOrig = ddq;
    
    q = [qOrig.ankle, qOrig.knee, qOrig.hip];
    dq = [dqOrig.ankle, dqOrig.knee, dqOrig.hip];
    ddq = [ddqOrig.ankle, ddqOrig.knee, ddqOrig.hip];
    
    tau = zeros(size(q));
    for indTime = 1:length(time.time) % recalc torque given redistributed masses
        %                 model.updateState(q(indTime, :), dq(indTime, :));
        tau(indTime, :) = model.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
    end
    
    trajT = time.time';
    trajX = encodeState(q, dq);
    trajU = tau;
    
    states = trajX;
    control = trajU;
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

function [progressVar] = calcWinLenAndH(trajT, trajX, trajU, dt, ioc, startInd, trajH1, trajH2, runParam)
    lenTraj = size(trajX, 1);
    
    H1 = [];
    H2 = [];
    prevFullWinInds = [];
    currFullWinInds = [];
    
    for i = 1:runParam.maxWinLen
        % determine the current window to check
        switch runParam.hWinAdvFlag
            case 'forward'
                currFullWinInds = startInd + (1:i);
                
            case 'backward'
                currFullWinInds = startInds - (1:i);
                
            case 'centre'
                addToRight = ceil((i-1)/2);
                addToLeft = ceil((i-2)/2);
                
                addInds = [-addToLeft:0 0:addToRight];
                uniInds = unique(addInds);
                currFullWinInds = uniInds + startInd;
        end
        
        % error check on the win len
        currFullWinInds = currFullWinInds(currFullWinInds > 0);
        currFullWinInds = currFullWinInds(currFullWinInds <= lenTraj);
        
        currLen = length(currFullWinInds);
        currHRow = currLen*runParam.lenDof;
        
        if currLen < runParam.winSize
            switch runParam.displayInfo
                case 'verbose'
                    fprintf("Obs frames %i:%i, width not sufficient. \n", ...
                        currFullWinInds(1), currFullWinInds(end));
            end
            continue;
            
        elseif currHRow < runParam.hWant % run the inversion only when the window is properly sized
            switch runParam.displayInfo
                case 'verbose'
                    fprintf("Obs frames %i:%i, row to col ratio not sufficient. Have %u but want %u \n", ...
                        currFullWinInds(1), currFullWinInds(end), currHRow, runParam.hWant);
            end
            continue;   
        end
        
        % check the previous entry, if the first and last entry matches,
        % then we can just reuse the previous H1/H2 matrices instead of
        % recon a new one
        if ~isempty(prevFullWinInds) && length(currFullWinInds) > 1 &&...
                prevFullWinInds(1) == currFullWinInds(1) && ...
                prevFullWinInds(end) == currFullWinInds(end-1)
            % use the previous H1 and H2
        else
            H1 = [];
            H2 = [];
        end
        
        % assemble H matrix
        [H, H1, H2] = assembleHMatrix(trajT, trajX, trajU, dt, ioc, H1, H2, currFullWinInds, trajH1, trajH2, runParam);
     
        prevFullWinInds = currFullWinInds;
        
        if currLen >= runParam.maxWinLen || max(currFullWinInds) == size(trajX, 1)
            hitMaxWinLen = 1;
        else
            hitMaxWinLen = 0;
        end
        
        % check H matrix for completion
        progressVar = checkHMatrix(H, currFullWinInds, runParam, hitMaxWinLen);
        
        if ~isempty(progressVar)
            H1 = [];
            H2 = [];
            break;
        end
    end
end

function [H, H1, H2] = assembleHMatrix(trajT, trajX, trajU, dt, ioc, H1, H2, fullWinInds, trajH1, trajH2, runParam)
    if ~isempty(H1) 
        indsToRun = length(fullWinInds);
    else
        indsToRun = 1:length(fullWinInds);
    end
    
    for indSubFrame = indsToRun
        % start each individual window
        winInds = fullWinInds(1:indSubFrame);

        % Read next observation
        x = trajX(winInds, :);
        u = trajU(winInds, :);

        % Assemble the H matrix for this width
        [H1, H2] = getRecoveryMatrix(ioc, H1, H2, x, u, dt);
        H = [H1 -H2];
    end
end

function progressVar = checkHMatrix(H, fullWinInds, runParam, hitMaxWinLen)
    progressVar = [];

    % proceed with H matrix inversion
    Hhat = H/norm(H,'fro');

    % Compute weights using final recovery matrix
    [weights, ~] = computeWeights(Hhat, runParam.numWeights);
    % weightTraj(indFrame, :) = weights/sum(weights);

    % Compute error between true and recovered weights
    error = computeError(weights, runParam.trueWeights);
%     errorTraj(indFrame,:) = error;

    % save rank of recovery matrix
    [rank, completed, rankPassCode] = validateCompletion(Hhat, runParam.gamma, runParam.delta, runParam.dimWeights);

    if hitMaxWinLen
        completed = 1;
    end

%     rankTraj(currLen,:) = rank;
%     rankPassCodeTraj(currLen, :) = rankPassCode;

switch runParam.displayInfo
    case 'verbose'
        % Display information
        fprintf("Obs frames %i:%i of max %u/%u, rank: %0.2f/%0.2f, error: %0.3f, code: %u,%u,%u ", ...
            fullWinInds(1), fullWinInds(end), length(fullWinInds), runParam.maxWinLen, rank, runParam.gamma, error(1), ...
            rankPassCode(1), rankPassCode(2), rankPassCode(3) );
        
        fprintf('Weights: ');
        for j = 1:runParam.numWeights
            fprintf('%+0.3f,', weights(j));
        end
        fprintf('\n');
        
    case 'final'
        if completed
            % Display information
            fprintf("Obs frames %i:%i of max %u/%u, rank: %0.2f/%0.2f, error: %0.3f, code: %u,%u,%u ", ...
                fullWinInds(1), fullWinInds(end), length(fullWinInds), runParam.maxWinLen, rank, runParam.gamma, error(1), ...
                rankPassCode(1), rankPassCode(2), rankPassCode(3) );
            
            fprintf('Weights: ');
            for j = 1:runParam.numWeights
                fprintf('%+0.3f,', weights(j));
            end
            fprintf('\n');
        end
end
        
    if completed
        progressVar.Hhat = Hhat;
        progressVar.weights = weights'/sum(weights);
        progressVar.winInds = fullWinInds;
        progressVar.rankTraj = rank;
        progressVar.rankPass = rankPassCode;
        progressVar.error = error;
    end
end