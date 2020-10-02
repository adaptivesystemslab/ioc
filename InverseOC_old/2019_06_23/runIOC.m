function ret = runIOC(trialInfo, savePath, displayInfo)

    ret = 0;

    % Check function arguments and define default values if needed
    if ~exist('displayInfo', 'var')
        displayInfo = 0;
    end
    
    if ~exist('savePath', 'var')
        savePath = "";
    end
    
    % Return to calling function if no trial info is passed along
    if isempty(fieldnames(trialInfo))
        disp("No trial info was given");
        return
    end
    
    % Load to memory model and information needed for running IOC
    model = getModel(trialInfo.model, trialInfo.modelType);
    trueWeights = (trialInfo.weights)'/sum(trialInfo.weights);
    numWeights = length(trueWeights);
    numDofs = model.totalActuatedJoints;
    
    % Create IOC instance
    ioc = IOCInstance(model, trialInfo.candidateFeatures);
    
    % Load state, control and time trajectories to be analyzed.
    inputPath = char(trialInfo.path);
    %inputPath = char(sprintf("%s%s", trialInfo.path, trialInfo.name));
    
    data = cell2mat(struct2cell(load(inputPath)));
    time = data(:,1);
    control = data(:,2:numDofs+1);
    stop = 10;
    states = data(:, numDofs+2:stop);

   % Approximate trajectories using splines
    trajT = linspace(0,1,1001);
    trajU = interp1(time, control, trajT,'spline');
    trajX = interp1(time, states, trajT,'spline');    
    
    % Initialization of threshold and condition variables
    start=2;    
    gamma = 100;
    delta=1e-6;
    dt = trajT(2) - trajT(1);
    dimWeights = 20;
    numObs = 0;
    
    % Initialization of concatenated arrays summarizing error and gamma as functions of t
    rankTraj = [];
    errorTraj = []; 
    
    % Build recovery matrix for first iteration
    x = trajX(start,:);
    u = trajU(start,:);
    x0 = trajX(start-1,:);
    u0 = trajU(start-1,:);

    [H1, H2] = getRecoveryMatrix(ioc, [], [], x0, u0, x, u, dt);
    H = [H1 -H2];
    Hhat = H/norm(H,'fro');
    
    for maxFrames = 1:size(trajT,2)-start
                                   
        % Move to next state
        numObs = numObs+1;
        next = start + numObs;
   
        % Read next observation
        x = trajX(next,:);
        u = trajU(next,:);
        x0 = trajX(next-1,:);
        u0 = trajU(next-1,:);
        
        % Generate new recovery matrix
        [H1, H2] = getRecoveryMatrix(ioc, H1, H2, x0, u0, x, u, dt);
        H = [H1 -H2];
        Hhat = H/norm(H,'fro');
    
        % Compute weights using final recovery matrix
        [weights, ~] = computeWeights(Hhat, numWeights);
        weights = weights/sum(weights);

        % Compute error between true and recovered weights
        error = computeError(weights, trueWeights);
        errorTraj(end+1,:) = error;
        
        % save rank of recovery matrix
        [rank, ~] = validateCompletion(Hhat, gamma, delta, dimWeights);
        rankTraj(end+1,:) = rank;

        if displayInfo
        % Display information
            fprintf("%i observations in recovery matrix, estimation error is %f\n", numObs, error(1));
        end
    end
    
    mat = [errorTraj, rankTraj];
    finalPath = sprintf("%s%s_%s.mat", savePath, trialInfo.name, trialInfo.model);
    save(char(finalPath), 'mat');
    ret = 1;

end

