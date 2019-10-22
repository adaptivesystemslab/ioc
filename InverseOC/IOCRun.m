function IOCRun(trialInfo, savePath)
    frameDelim = '_';
    
    % Load to memory model and information needed for running IOC
    model = getModel(trialInfo);
    traj = model.loadData(trialInfo);
    
    q = traj.q;
    dq = traj.dq;
    tau = traj.tau;
    frameInds = traj.frameInds;
    trajT = traj.trajT;
    trajX = traj.trajX;
    trajU = traj.trajU;
    
    if 0 % visualize the data if want to check the q generated
        model.plotTrajectory(q);
    end

    trialInfo.frameInds = frameInds;
    trialInfo.numWeights = length(trialInfo.candidateFeatures);
    trialInfo.numDofs = model.getModelDof;
    
    % generate default trialInfo if missing weights
    if isempty(trialInfo.weights)
        trialInfo.weights = ones(1, trialInfo.numWeights);
    end
    
    % gamma
    if isempty(trialInfo.gamma)
        trialInfo.gamma = trialInfo.numWeights+trialInfo.numDofs-1;
    end
    
    trialInfo.trueWeights = (trialInfo.weights)'/sum(trialInfo.weights);
    
    % Initialization of threshold and condition variables
    dt = trajT(2) - trajT(1);

    % Create IOC instance
    ioc = IOCInstance(model, dt);
    ioc.init(trialInfo);
%     ioc.setFeatureNormalization(trajX, trajU);
    
    switch trialInfo.displayInfo
        case {'verbose', 'final'}
            fprintf('Pre-calculating the features/dynamics... \n');
    end
        
    precalcAllFrames = [frameInds frameInds(end)+(1:trialInfo.maxWinLen)];
    
    if max(precalcAllFrames) > size(trajX, 1)
        precalcAllFrames = precalcAllFrames(1):size(trajX, 1);
    end
    
    trialInfo.lenDof = size(trajX, 2)/2;
    trialInfo.hWant = (size(trajX, 2) + trialInfo.numWeights) * trialInfo.dimWeights;
    
    % calculating the features/dynamics for storage prurposes
    iocFeatures = ioc.calcFeatures(trajX(precalcAllFrames, :), trajU(precalcAllFrames, :));
    iocDynamics = ioc.calcDynamics(trajX(precalcAllFrames, :), trajU(precalcAllFrames, :));
 
    if trialInfo.saveIntermediate > 0 % save a mat file based on the value here, to limit the size of the overall mat file
        saveInds = [];
        for i = 1:ceil(frameInds(end)/trialInfo.saveIntermediate)
            saveInds = [saveInds; (i-1)*trialInfo.saveIntermediate i*trialInfo.saveIntermediate-1];
        end
        
        saveInds(1, 1) = 1;
        saveInds(end, 2) = frameInds(end);
    else
        saveInds(1, 1) = 1;
        saveInds(1, 2) = frameInds(end);
    end
    
    % print to console so it will be saved in output logs
    trialInfo
    
    % precalc singular H1 and H2 matrix
    precalcGradient = precalculateGradient_initialize(trajX, trajU, ioc, 1:trialInfo.maxWinLen, trialInfo);
    
    progressVar = [];
    progressVar(frameInds(end)).weights = [];
    progressVar(frameInds(end)).winInds = [];
    progressVar(frameInds(end)).rankTraj = [];
    progressVar(frameInds(end)).rankPass = [];
    progressVar(frameInds(end)).error = [];
    
    processSecondaryVar = [];
    processSecondaryVar(frameInds(end)).H1 = [];
    processSecondaryVar(frameInds(end)).H2 = [];
    processSecondaryVar(frameInds(end)).H = [];
    
    % assemble and assess H
    for indFrame = frameInds
        [progressVarTemp, processSecondaryVarTemp, precalcGradient] = calcWinLenAndH(trajT, trajX, trajU, dt, ioc, indFrame - 1, precalcGradient, trialInfo);
        
        if ~isempty(progressVarTemp)
            progressVar(indFrame) = progressVarTemp;
            processSecondaryVar(indFrame) = processSecondaryVarTemp;
        end
        
        % save either on the regular intervals, or at the last frame
        checkSave = find(indFrame == saveInds(:, 2));
        
        if (~isempty(checkSave))    
            currSaveInds = saveInds(checkSave, :);
            currSaveRange = currSaveInds(1):currSaveInds(end);
            
            outputVar_data.frameInds = currSaveInds;
            outputVar_data.t = trajT(currSaveRange);
            outputVar_data.q = q(currSaveRange, :);
            outputVar_data.dq = dq(currSaveRange, :);
            outputVar_data.tau = tau(currSaveRange, :);
            outputVar_data.dt = dt;
            
            outputVar_data.lenDof = trialInfo.lenDof;
            outputVar_data.lenState = size(trajX, 2);
            outputVar_data.lenControl = size(trajU, 2);
            
            outputVar_data.featureLabels = ioc.getFeatureListAsString();
            outputVar_data.features = iocFeatures(currSaveRange, :);
            outputVar_data.dynamics = iocDynamics(currSaveRange, :);
                        
            outputVar_weights.trialInfo = trialInfo;
            outputVar_weights.timeElapsed = toc;
            outputVar_weights.timeElasedPerFrame = toc/length(frameInds);
            
            outputVar_weights.frameInds = currSaveInds;
            outputVar_weights.progress = progressVar(currSaveRange);
            
            outputVar_weights.minLenThres = trialInfo.hWant / outputVar_data.lenDof;
            outputVar_weights.maxLenThres = trialInfo.maxWinLen;
            outputVar_weights.minRankThres = trialInfo.gamma;
            
            % outputVar.errorTraj = errorTraj;
            % outputVar.rankTraj = rankTraj;
            % outputVar.weightTraj = weightTraj;
            % outputVar.completeTraj = completeTraj;
            % outputVar.rankPassCodeTraj = rankPassCodeTraj;
            
            % outputVar.segmentArray = segmentArray; % remove the first initialized value
            outputVar_supp.frameInds = currSaveInds;
            outputVar_supp.processSecondaryVar = processSecondaryVar(currSaveRange);
            
            numSuffix = [num2str(currSaveInds(1), '%06.f') '_' num2str(currSaveInds(2), '%06.f')];

            finalPath_data = fullfile(savePath, ['data' frameDelim numSuffix '.mat']);
            finalPath_weights = fullfile(savePath, ['weights' frameDelim numSuffix '.mat']);
            finalPath_supp = fullfile(savePath, ['supp' frameDelim numSuffix '.mat']); 
            
            save(finalPath_data, 'outputVar_data');
            save(finalPath_weights, 'outputVar_weights');
            save(finalPath_supp, 'outputVar_supp');
            
            progressVar = [];
            progressVar(frameInds(end)).weights = [];
            progressVar(frameInds(end)).winInds = [];
            progressVar(frameInds(end)).rankTraj = [];
            progressVar(frameInds(end)).rankPass = [];
            progressVar(frameInds(end)).error = [];
            
            processSecondaryVar = [];
            processSecondaryVar(frameInds(end)).H1 = [];
            processSecondaryVar(frameInds(end)).H2 = [];
            processSecondaryVar(frameInds(end)).H = [];
        end
    end
end

function [progressVar, processSecondaryVar, precalcGradient] = calcWinLenAndH(trajT, trajX, trajU, dt, ioc, startInd, precalcGradient, trialInfo)
    lenTraj = size(trajX, 1);
    
    progressVar = [];
    processSecondaryVar = [];
    
    H1 = [];
    H2 = [];
    prevFullWinInds = [];
    currFullWinInds = [];
    
    % frameIndsFullRange
    frameIndsFullRange = (startInd-trialInfo.maxWinLen):(startInd+trialInfo.maxWinLen);
    frameIndsFullRange = frameIndsFullRange(frameIndsFullRange > 0);
    frameIndsFullRange = frameIndsFullRange(frameIndsFullRange <= trialInfo.frameInds(end)+trialInfo.maxWinLen);
    frameIndsFullRange = frameIndsFullRange(frameIndsFullRange <= length(trajT));
    
    precalcGradient = precalculateGradient_pushpop(trajX, trajU, ioc, precalcGradient, frameIndsFullRange);
    
    for i = 1:trialInfo.maxWinLen
        % determine the current window to check
        switch trialInfo.hWinAdvFlag
            case 'forward'
                currFullWinInds = startInd + (1:i);
                
            case 'backward'
                currFullWinInds = startInd - (1:i);
                
            case 'centre'
                addToRight = ceil((i-1)/2);
                addToLeft = ceil((i-2)/2);
                
                addInds = [-addToLeft:0 0:addToRight];
                uniInds = unique(addInds);
                currFullWinInds = uniInds + startInd + 1;
        end
        
        % error check on the win len
        currFullWinInds = currFullWinInds(currFullWinInds > 0);
        currFullWinInds = currFullWinInds(currFullWinInds <= lenTraj);
        
        currLen = length(currFullWinInds);
        currHRow = currLen*trialInfo.lenDof;
        
        if currLen < trialInfo.windowWidth
            switch trialInfo.displayInfo
                case 'verbose'
                    fprintf('Obs frames %i:%i, width not sufficient. \n', ...
                        currFullWinInds(1), currFullWinInds(end));
            end
            continue;
            
        elseif currHRow < trialInfo.hWant % run the inversion only when the window is properly sized
            switch trialInfo.displayInfo
                case 'verbose'
                    fprintf('Obs frames %i:%i, row to col ratio not sufficient. Have %u but want %u \n', ...
                        currFullWinInds(1), currFullWinInds(end), currHRow, trialInfo.hWant);
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
        [H, H1, H2] = assembleHMatrixWithPrecalc(H1, H2, currFullWinInds, precalcGradient);
%         [H, H1, H2] = assembleHMatrix(trajT, trajX, trajU, dt, ioc, H1, H2, currFullWinInds, trajH1, trajH2, trialInfo);
     
        prevFullWinInds = currFullWinInds;
        
        if currLen >= trialInfo.maxWinLen || max(currFullWinInds) == size(trajX, 1)
            hitMaxWinLen = 1;
        else
            hitMaxWinLen = 0;
        end
        
        % check H matrix for completion
        [progressVar] = checkHMatrix(H, currFullWinInds, trialInfo, trajT, hitMaxWinLen);
        
        if ~isempty(progressVar)
            processSecondaryVar.H1 = H1;
            processSecondaryVar.H2 = H2;
            processSecondaryVar.H = H;
%             processSecondaryVar.Hhat = Hhat; %  Hhat = H/norm(H,'fro');
            
            break;
        end
    end
end

function [H, H1, H2] = assembleHMatrix(trajT, trajX, trajU, dt, ioc, H1, H2, fullWinInds, trialInfo)
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

function [df_dx, df_du, dp_dx, dp_du, x, u] = getGradient(trajX, trajU, incurrInd, ioc)   
    [x, u] = setDataLength(trajX, trajU, IOCInstanceNew.winSize, incurrInd);
    [fx, fu, px, pu] = ioc.getDerivativesNewObservation(x, u);
    df_dx = fx';
    df_du = fu';
    dp_dx = px';
    dp_du = pu'; 
end

function [H, H1, H2] = assembleHMatrixWithPrecalc(H1, H2, fullWinInds, precalcGradient)
    % need to recalculate H1/H2 using precalcs to the current timestep
    if isempty(H1)
        for indSubFrame = 1:length(fullWinInds)-1
            currFrame = fullWinInds(indSubFrame);
            [H1, H2] = assembleH1H2(precalcGradient(currFrame).df_dx, ...
                precalcGradient(currFrame).df_du, ...
                precalcGradient(currFrame).dp_dx, ...
                precalcGradient(currFrame).dp_du, H1, H2);
        end
    end
    
    % add the newest timestep to it
    currFrame = fullWinInds(end);
    [H1, H2] = assembleH1H2(precalcGradient(currFrame).df_dx, ...
        precalcGradient(currFrame).df_du, ...
        precalcGradient(currFrame).dp_dx, ...
        precalcGradient(currFrame).dp_du, H1, H2);

    H = [H1 -H2];
end

function [progressVar] = checkHMatrix(H, fullWinInds, trialInfo, trajT, hitMaxWinLen)
    progressVar = [];

    % proceed with H matrix inversion
    Hhat = H/norm(H,'fro');

    % Compute weights using final recovery matrix
    [weights, ~] = computeWeights(Hhat, trialInfo.numWeights);
    % weightTraj(indFrame, :) = weights/sum(weights);

    % Compute error between true and recovered weights
    error = computeError(weights, trialInfo.trueWeights);
%     errorTraj(indFrame,:) = error;

    % save rank of recovery matrix
    [rank, completed, rankPassCode] = validateCompletion(Hhat, trialInfo.gamma, trialInfo.delta, trialInfo.dimWeights);

    if hitMaxWinLen
        completed = 1;
    end

%     rankTraj(currLen,:) = rank;
%     rankPassCodeTraj(currLen, :) = rankPassCode;

    switch trialInfo.displayInfo
        case 'verbose'
            % Display information
% % %             fprintf('Obs frames %i:%i of max %u/%u, rank: %0.2f/%0.2f, error: %0.3f, code: %u,%u,%u ', ...
% % %                 fullWinInds(1), fullWinInds(end), length(fullWinInds), trialInfo.maxWinLen, rank, trialInfo.gamma, error(1), ...
% % %                 rankPassCode(1), rankPassCode(2), rankPassCode(3) );

           fprintf('[%s] Obs %i:%i of %i, len %u/%u, rank %0.2f/%0.2f, ', ...
                datestr(now), fullWinInds(1), fullWinInds(end), length(trajT), length(fullWinInds), trialInfo.maxWinLen, rank, trialInfo.gamma);

            fprintf('Weights: ');
            for j = 1:trialInfo.numWeights
                fprintf('%+0.3f,', weights(j));
            end
            fprintf('\n');

        case 'final'
            if completed
                % Display information
                fprintf('[%s] Obs %i:%i of %i, len %u/%u, rank %0.2f/%0.2f, ', ...
                    datestr(now), fullWinInds(1), fullWinInds(end), length(trajT), length(fullWinInds), trialInfo.maxWinLen, rank, trialInfo.gamma);
                
                fprintf('Weights: ');
                for j = 1:trialInfo.numWeights
                    fprintf('%+0.3f,', weights(j));
                end
                fprintf('\n');
            end
    end
        
    if completed
%         progressVar.Hhat = Hhat;
        progressVar.weights = weights'/sum(weights);
%         progressVar.winIndsStart = fullWinInds(1);
%         progressVar.winIndsEnd = fullWinInds(end);
        progressVar.winInds = [fullWinInds(1) fullWinInds(end)];
        progressVar.rankTraj = rank;
        progressVar.rankPass = rankPassCode;
        progressVar.error = error;
    end
end

function precalcGradient = precalculateGradient_initialize(trajX, trajU, ioc, frameInds, trialInfo)
    precalcGradient(frameInds(end)+trialInfo.maxWinLen).df_dx = []; % preallocating
    precalcGradient(frameInds(end)+trialInfo.maxWinLen).df_du = [];
    precalcGradient(frameInds(end)+trialInfo.maxWinLen).dp_dx = [];
    precalcGradient(frameInds(end)+trialInfo.maxWinLen).dp_du = [];
    
    for i = frameInds
        fprintf('Pre-calculating %uth H1/H2... \n', i);
        [tempGrad.df_dx, tempGrad.df_du, tempGrad.dp_dx, tempGrad.dp_du] = ...
            getGradient(trajX, trajU, i, ioc);
        
        precalcGradient(i) = tempGrad;
    end
end

function precalcGradient = precalculateGradient_pushpop(trajX, trajU, ioc, precalcGradient, frameInds)
    % pop the frame before
    priorInd = frameInds(1) - 1;
    if priorInd > 0
        precalcGradient(priorInd).df_dx = []; % preallocating
        precalcGradient(priorInd).df_du = [];
        precalcGradient(priorInd).dp_dx = [];
        precalcGradient(priorInd).dp_du = [];
    end
    
    % generate the new frame
    newFrameInd = frameInds(end)+1;
   
    if newFrameInd <= length(trajX)
        [tempGrad.df_dx, tempGrad.df_du, tempGrad.dp_dx, tempGrad.dp_du] = ...
            getGradient(trajX, trajU, newFrameInd, ioc);
        
        precalcGradient(newFrameInd) = tempGrad;
        
% % %         fprintf('Removed %uth H1/H2 and adding %uth H1/H2... ', priorInd, newFrameInd);
    end
end