function matData = assembleData(currSubPath)
    dirInnerPath = dir(fullfile(currSubPath, 'weights*.mat'));
    
    if isempty(dirInnerPath)
        matData = [];
        return;
    end
    
    for i = 1:length(dirInnerPath)
        currFileName = dirInnerPath(i).name;
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        currInds = outputVar_weights.frameInds(1):outputVar_weights.frameInds(end);
        matData.timeElapsed(i) = outputVar_weights.timeElapsed;
        matData.progress(currInds) = outputVar_weights.progress;
        
        currFileName = strrep(dirInnerPath(i).name, 'weights', 'data');
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        matData.t(currInds) = outputVar_data.t;
        matData.q(currInds, :) = outputVar_data.q;
%         matData.dq(currInds, :) = outputVar_data.dq;
%         matData.tau(currInds, :) = outputVar_data.tau;
%         matData.features(currInds, :) = outputVar_data.features;
%         matData.dynamics(currInds, :) = outputVar_data.dynamics;

        currFileName = strrep(dirInnerPath(i).name, 'weights', 'supp');
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        matData.processSecondaryVar(currInds) = outputVar_supp.processSecondaryVar;
    end
    
    matData.trialInfo = outputVar_weights.trialInfo;
    matData.timeElapsed = max(matData.timeElapsed);
    matData.minLenThres = outputVar_weights.minLenThres;
    matData.maxLenThres = outputVar_weights.maxLenThres;
    matData.minRankThres = outputVar_weights.minRankThres;
    
    matData.lenDof = outputVar_data.lenDof;
    matData.lenState = outputVar_data.lenState;
    matData.lenControl = outputVar_data.lenControl;
    matData.featureLabels = outputVar_data.featureLabels;
    matData.dt = outputVar_data.dt;
    
    % now calculate a few things
    matData.lenWeights = length(matData.featureLabels);
end