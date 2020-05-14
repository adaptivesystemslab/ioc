function templateModel = setupForFg(featureModel)
    templateModel.featureModel = featureModel;
    
    % create empty HMMs since we don't have HMM models
    hmmModel.LLThresholdAbs = 0;
    hmmModel.LLThreshold = 0;
    
    for i = 1:length(featureModel)
        templateModel.hmmModel{i} = hmmModel;
        templateModel.motionName{i} = featureModel{i}.name;
        
        if length(featureModel{i}.veloTime{1}) > 1
            meanVeloTime(i) = diff(featureModel{i}.veloTime{1});
        else
            meanVeloTime(i) = featureModel{i}.veloTime{1};
        end
    end
    
    templateModel.meanPeakToPeak = mean(meanVeloTime);
    
    templateModel.minTemplateTime = 3;
    templateModel.allSigDof = 1;
    templateModel.lenTemplate = length(featureModel);