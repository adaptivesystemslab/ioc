function ret = runIOC(trialInfo, savePath, runParam)
    % Check function arguments and define default values if needed
    displayInfo = runParam.displayInfo;
    saveIntermediate = runParam.saveIntermediate;
    gamma = runParam.gamma;
    dimWeights = runParam.dimWeights;
    start = runParam.start;
    
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
    numWeights = length(trialInfo.candidateFeatures);
    numDofs = model.totalActuatedJoints;
    trueWeights = (trialInfo.weights)'/sum(trialInfo.weights);
    
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
     
    % Initialization of threshold and condition variables
    numObs = 0;
    completed = 1;
    dt = trajT(2) - trajT(1);

    % Create IOC instance
    ioc = IOCInstanceNew(model, dt);
    ioc.init(trialInfo);
    ioc.setFeatureNormalization(trajX, trajU);
    
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
        mdl.forwardPosition();
        vis.update();
        
        for s = 1:length(trajT)
            mdl.position = q(s, :);
            mdl.forwardPosition();
            vis.update();
        end
    end

    % Initialization of concatenated arrays summarizing error and gamma as functions of t
    frameInds = 1:size(trajT,2);
    weightTraj = zeros(size(frameInds, 2), numWeights);
    errorTraj = zeros(size(frameInds))';
    rankTraj = zeros(size(frameInds))';
    completeTraj = zeros(size(frameInds))';
    rankPassCodeTraj = zeros(length(frameInds), 3);
    segmentArray = []; 
    segmentArrayInds = 0;
    H1 = [];
    H2 = [];
    
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
        end
        
        % Move to next state
        numObs = numObs+1;
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
        [H1, H2] = getRecoveryMatrix(ioc, H1, H2, x0, u0, x, u, dt);
        H = [H1 -H2];
        
        [checkHbool, checkHhave, checkHwant] = hPassFct1(H, dimWeights);
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
            rankTraj(indFrame,:) = rank;
            completeTraj(indFrame, :) = completed;
            rankPassCodeTraj(indFrame, :) = rankPassCode;
            
            if displayInfo
                % Display information
                fprintf("(%uth seg) Obs len: %i/%i, rank: %0.2f/%0.2f, error: %0.3f, code: %u,%u,%u ", ...
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
            outputVar.t = time;
            outputVar.q = q;
            outputVar.dq = dq;
            outputVar.tau = tau;
            
            outputVar.featureLabels = ioc.getFeatureListAsString();
            outputVar.features = iocFeatures;
            outputVar.dynamics = iocDynamics;
            
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
