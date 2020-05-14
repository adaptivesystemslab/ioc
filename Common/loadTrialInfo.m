function [trialInfo] = loadTrialInfo(trialInfoOrig, configFile, potentialBasePaths, configFilePath, templateFile)
    for j = 1:length(potentialBasePaths)
        targetPath = fullfile(potentialBasePaths{j}, trialInfoOrig.subpath);
        
        if exist(targetPath, 'file')
            break;
        end
    end
    
    trialInfo = struct;
    % pre-populate the trialInfo with the global settings in the json file
    runParamFields = fieldnames(configFile.runParamGlobal);
    for j = 1:length(runParamFields)
        trialInfo.(runParamFields{j}) = configFile.runParamGlobal.(runParamFields{j});
    end
    
    % then search for the specific template settings that we're
    % interested. found, overwrite any globalsettings with the template
    for j = 1:length(templateFile.runTemplates)
        if strcmpi(templateFile.runTemplates(j).templateName, trialInfoOrig.runTemplate)
            runParam = templateFile.runTemplates(j);
            runParamFields = fieldnames(runParam);
            for k = 1:length(runParamFields)
                trialInfo.(runParamFields{k}) = runParam.(runParamFields{k});
            end
            break;
        end
    end

    % then overwrite any globalsettings with the specific settings in the
    % last block of the json file
    runParamFields = fieldnames(trialInfoOrig);    
    for j = 1:length(runParamFields) 
        trialInfo.(runParamFields{j}) = trialInfoOrig.(runParamFields{j});
    end
    
    trialInfo.configFile = configFilePath;
    trialInfo.windowWidth = IOCInstance.winSize;
    trialInfo.delta = IOCInstance.delta;
    trialInfo.path = targetPath;
end