function ret = runIOC_WJTest(trialInfo, savePath, runParam)
    % Check function arguments and define default values if needed
    displayInfo = runParam.displayInfo;
    saveIntermediate = runParam.saveIntermediate;
    gamma = runParam.gamma;
    maxWinLen = runParam.maxWinLen;
    dimWeights = runParam.dimWeights;
    start = runParam.start;
    
    winSize = IOCInstanceNew.winSize;
    delta = IOCInstanceNew.delta;
    
    frameInds = 1:1000;

    
    if ~exist('savePath', 'var')
        savePath = "";
    end
    
    % Return to calling function if no trial info is passed along
    if isempty(fieldnames(trialInfo))
        disp("No trial info was given");
        return
    end
    
    % Load to memory model and information needed for running IOC
%     model = getModel(trialInfo.model, trialInfo.modelType, [], trialInfo.path);
    model = ArmModelRL();
    model.loadModelOnly();
%     model.model.base = 'frame_6dof_root';
    model.forwardKinematics();
    
    numWeights = length(trialInfo.candidateFeatures);
    numDofs = model.totalActuatedJoints;
    trueWeights = (trialInfo.weights)'/sum(trialInfo.weights);
    
%     load('D:\aslab\projects\jf2lin\TROcopy\data\original\SQUA_STD_FAT_Subject01_export.mat');
%     ts=500;
%     te=7637;
%     q(:, 1) = saveVar.jointAngle.joint_rhip_0(ts:te);
%     q(:, 2) = saveVar.jointAngle.joint_rknee_0(ts:te);
%     q(:, 3) = saveVar.jointAngle.joint_rankle_0(ts:te)+pi/2;
    
    
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
    
    mdl = model.model;
%     mdl.g = [0 0 9.81];
    
    if 0
        qTest = q(1, :)
%         dqTest = dq(1, :);
%         ddqTest = ddq(1, :);

%         qTest = [0 0 -pi/2]
%         qTest = [0 0 0]
%         dqTest = zeros(3, 1);
%         ddqTest = zeros(3, 1);

        qTest = [1 1 1];
        dqTest = [2 2 2];
        ddqTest = [3 3 3];

        mdl.position = qTest;
        mdl.velocity = dqTest;
        mdl.acceleration = ddqTest;
        mdl.forwardPosition();mdl.forwardVelocity();mdl.forwardAcceleration();
        mdl.inverseDynamics();
        tau = mdl.torque'
        
        mdl.calculateMassMatrix();
        M = mdl.M;
        mdl.calculateGravity();
        G = mdl.G;
        
        mdl.position = qTest;
        mdl.velocity = dqTest;
        mdl.acceleration = ddqTest;
        mdl.forwardPosition();mdl.forwardVelocity();mdl.forwardAcceleration();
        mdl.inverseDynamics();
        mdl.calculateCentrifugalCoriolis();
        C = mdl.V;
        
        clc
        M*ddqTest'
        C
        G
        
        output = [M*ddqTest' C G]
    end

    tauRaw = zeros(size(q));
    for indTime = 1:length(time.time) % recalc torque given redistributed masses
        %                 model.updateState(q(indTime, :), dq(indTime, :));
        tauRaw(indTime, :) = model.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
    end
    
    
    tau = tauRaw; % [T.ankle, T.knee, T.hip];
    
    trajT = time.time';
    trajX = encodeState(q, dq);
    trajU = tau;

    % Initialization of threshold and condition variables
    numObs = 0;
    completed = 1;
    dt = trajT(2) - trajT(1);

    % Create IOC instance
    ioc = IOCInstanceNew(model, dt);
    ioc.init(trialInfo);
%     ioc.setFeatureNormalization(trajX, trajU);
    
    iocFeatures = ioc.calcFeatures(trajX, trajU);
    iocDynamics = ioc.calcDynamics(trajX, trajU);

    if 0
        % check the ddq/torque are internally consistent
        %          for indTime = 1:length(time)
        %              model.updateState(q(indTime, :), dq(indTime, :));
        %              ddqCheck(indTime, :) = model.forwardDynamics(tauRaw(indTime, :));
        %          end
        %
        % plot the q/dq/ddq/tau
%         figure;
%         ax(1) = subplot(221); plot(q); title('q');
%         ax(2) = subplot(222); plot(dq); title('dq');
%         ax(3) = subplot(223); plot(ddq); title('ddq');
%         ax(4) = subplot(224); plot(tau); title('tau');
%         linkaxes(ax, 'x');
        
%         % check the feature/dynamics output
%         f = ioc.calcFeatures(states, control);
%         p = ioc.calcDynamics(states, control);
%         [df_dx, df_du, dp_dx, dp_du] = ioc.getDerivativesNewObservation(trajX, trajU);

        featureLabels = ioc.getFeatureListAsString();
        subplot(2, 4, 1);
        plot(q); title('q');
        subplot(2, 4, 2);
        plot(dq); title('dq');
        subplot(2, 4, 3);
        plot(tau); title('tau');
        for i = 1:length(featureLabels)
            subplot(2, 4, 3+i);
            plot(iocFeatures(:, i));
            title(featureLabels{i});
        end
        
        mdl = model.model;
        vis = rlVisualizer('vis',640,480);
        mdl.forwardPosition();
        vis.addModel(mdl);
        vis.update();
        
        mdl.position = q(1, :);
%         mdl.position = [pi/2 0 0]
        mdl.forwardPosition();
        mdl.forwardVelocity();
        mdl.forwardAcceleration();
        mdl.inverseDynamics();
        vis.addMarker('x-axis', [1 0 0], [0.2 0.2 0.2 1]);
        vis.addMarker('y-axis', [0 1 0], [0.2 0.2 0.2 1]);
        vis.addMarker('z-axis', [0 0 1], [0.2 0.2 0.2 1]);
        for i = 1:length(mdl.joints)
            vis.addMarker(['j_' mdl.joints(i).name], mdl.joints(i).t(1:3, 4), [0.4 0.4 0.8 1]);
        end
        for i = 1:length(mdl.bodies)
            vis.addMarker(['EF_' mdl.bodies(i).name], mdl.bodies(i).t(1:3, 4), [0.4 0.8 0 1]);
            vis.addMarker(['COM_' mdl.bodies(i).name], mdl.bodies(i).t(1:3, 4) + mdl.bodies(i).com, [0.8 0 0.4 1]);
        end
        vis.update();
        
        for s = 1:8000
            mdl.position = q(s, :);
            mdl.forwardPosition();
            
            for i = 1:length(mdl.joints)
                vis.addMarker(['j_' mdl.joints(i).name], mdl.joints(i).t(1:3, 4), [0.4 0.4 0.8 1]);
            end
            for i = 1:length(mdl.bodies)
                vis.addMarker(['EF_' mdl.bodies(i).name], mdl.bodies(i).t(1:3, 4), [0.4 0.8 0 1]);
                vis.addMarker(['COM_' mdl.bodies(i).name], mdl.bodies(i).t(1:3, 4) + mdl.bodies(i).com, [0.8 0 0.4 1]);
            end
            
            vis.update();
            pause(0.01);
        end
    end

    % Initialization of concatenated arrays summarizing error and gamma as functions of t
    if isempty(frameInds)
        frameInds = 1:size(trajT,2);
    end
    
    weightTraj = zeros(size(frameInds, 2), numWeights);
    errorTraj = zeros(size(frameInds))';
    rankTraj = zeros(size(frameInds))';
    completeTraj = zeros(size(frameInds))';
    rankPassCodeTraj = zeros(length(frameInds), 3);
    segmentArray = []; 
    segmentArrayInds = 0;
    H1 = [];
    H2 = [];
    currLen = 0;
    
    for indFrame = frameInds
        if completed
            segmentArrayInds = segmentArrayInds + 1;          
            if segmentArrayInds > 1
                % if there are existing segment array, populate the
                % previous time step
                segmentArray(segmentArrayInds-1, 2) = indFrame - 1;
            end
            
            % populate the new one
            segmentArray(segmentArrayInds, 1) = indFrame;
            
            H1Log{segmentArrayInds} = H1;
            H2Log{segmentArrayInds} = H2;
            
            H1 = [];
            H2 = [];
            completed = 0;
            currLen = 0;
        end
        
        % Move to next state
        numObs = numObs+1;
        currLen = currLen+1;
        nextInds = (-winSize:0)+(start+numObs);
        next0Inds = (-winSize:0)+(start+numObs-1);
        
        if min(next0Inds) < 1
            if displayInfo
                fprintf("Obs len: %i/%i, width not sufficient. Min index is %u \n", ...
                    numObs, length(frameInds), min(next0Inds));
            end
            continue;
        end
        
        % Read next observation
        x = trajX(nextInds,:);
        u = trajU(nextInds,:);
        x0 = trajX(next0Inds,:);
        u0 = trajU(next0Inds,:);
        
        % Generate new recovery matrix
        [H1, H2] = getRecoveryMatrix(ioc, H1, H2, x0, u0, dt);
        H = [H1 -H2];
        
%         [checkHbool, checkHhave, checkHwant] = hPassFct1(H, dimWeights); 
        checkHbool = 1;
        checkHhave = 1;
        checkHwant = 1;
        if ~checkHbool % run the inversion only when the window is properly sized
            if displayInfo
                fprintf("Obs len: %i/%i, row to col ratio not sufficient. Have %u but want %u \n", ...
                    numObs, length(frameInds), checkHhave, checkHwant);
            end
        else
            Hhat = H/norm(H,'fro');
            
            % Compute weights using final recovery matrix
            [weights, ~] = computeWeights(Hhat, numWeights);
            weightTraj(indFrame, :) = weights/sum(weights);
            
            % Compute error between true and recovered weights
            error = computeError(weights, trueWeights);
            errorTraj(indFrame,:) = error;
            
            % save rank of recovery matrix
            [rank, completed, rankPassCode] = validateCompletion(Hhat, gamma, delta, dimWeights);
            
            if currLen > maxWinLen
                completed = 1;
            end
            
            rankTraj(indFrame,:) = rank;
            completeTraj(indFrame, :) = completed;
            rankPassCodeTraj(indFrame, :) = rankPassCode;
            
            if displayInfo
                % Display information
                fprintf("(Seg %u) Obs len: %i/%i, rank: %0.2f/%0.2f, error: %0.3f, code: %u,%u,%u ", ...
                    segmentArrayInds, numObs, length(frameInds), rank, gamma, error(1), rankPassCode(1), rankPassCode(2), rankPassCode(3) );
                
                fprintf('Weights: ');
                for j = 1:numWeights
                    fprintf('%+0.3f,', weights(j));
                end
                fprintf('\n');
            end
        end
        
        % save either on the regular intervals, or at the last frame
        if (saveIntermediate > 0 && mod(indFrame, saveIntermediate) == 0) || indFrame == frameInds(end)
            outputVar.t = trajT;
            outputVar.q = q;
            outputVar.dq = dq;
            outputVar.tau = tau;
            
            outputVar.featureLabels = ioc.getFeatureListAsString();
            outputVar.features = iocFeatures;
            outputVar.dynamics = iocDynamics;
            
            outputVar.frameInds = frameInds;
            
            outputVar.runParam = runParam;
            outputVar.trialInfo = trialInfo;
            
            outputVar.errorTraj = errorTraj;
            outputVar.rankTraj = rankTraj;
            outputVar.weightTraj = weightTraj;
            outputVar.completeTraj = completeTraj;
            outputVar.rankPassCodeTraj = rankPassCodeTraj;
            outputVar.timeElapsed = toc;
            outputVar.segmentArray = segmentArray; % remove the first initialized value
            outputVar.H1Log = H1Log;
            outputVar.H2Log = H2Log;
            
            finalPath = sprintf("%s%s_%s_struct.mat", savePath, trialInfo.name, trialInfo.model);
            save(char(finalPath), 'outputVar');
        end
    end
end
