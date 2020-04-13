function IOCAnalysis()
    setPaths();
    
    faceColours = brewermap(8, 'paired');
    
    nowstr = datestr(now, 'yyyymmddHHMMSS');
%     inputBasePath = 'D:\aslab_github\ioc\Data\IOC\';
%     outputBasePath = ['D:\aslab_github\ioc\Data\IOC\plot_' nowstr];
    inputBasePath = 'D:\aslab_github\ioc\Data\IOC\';
    outputBasePath = ['D:\aslab_github\ioc\Data\IOC\demo_plot\' nowstr];
    masterPathCsv = [outputBasePath '\a_summary.csv'];
    checkMkdir(outputBasePath);
    
    % iterate through all the folders in bathpath to produce figures and
    % table summary of the results
    loadPathInd = 0;
    dirPath = dir(inputBasePath);
    for i = 1:length(dirPath)
        currPath = fullfile(inputBasePath, dirPath(i).name);
        if strcmpi(currPath(end), '.')
            continue; % it's . or ..
        end
        
        dirCurrSubPath = dir(currPath);
        for j = 1:length(dirCurrSubPath)
            currSubPath = fullfile(currPath, dirCurrSubPath(j).name);
            if strcmpi(currSubPath(end), '.')
                continue; % it's . or ..
            end
            
            if ~exist(currSubPath, 'dir')
                continue;
            end
           
            loadPathInd = loadPathInd + 1;
            loadPaths{loadPathInd} = currSubPath;
        end
    end
    
    for i = 1:length(loadPaths)
        matData = assembleData(currSubPath);
        
        if isempty(matData)
            continue;
        end
        
        suffix = [dirPath(i).name '_' matData.trialInfo.runName '_' matData.trialInfo.templateName];
        
        outputPath = [outputBasePath '\' dirPath(i).name];
        checkMkdir(outputPath);
        
        try
            % plot results
            fprintf('%s\n', suffix);
            outputPathFig1 = fullfile(outputPath, ['fig_results_individual_' suffix]);
            outputPathFig2 = fullfile(outputPath, ['fig_results_cumulativeAllPass_' suffix]);
            outputPathFig3 = fullfile(outputPath, ['fig_results_cumulativeRankPass_' suffix]);
            outputPathMat1 = fullfile(outputPath, ['mat_results_cumulativeAllPass_' suffix]);
            outputPathCsv = fullfile(outputPath, ['csv_' suffix]);
            
            plotting_individual(matData, outputPathFig1, outputPathCsv, masterPathCsv);
            matSave = plotting_cumulative(matData, outputPathFig2, outputPathFig3, outputPathCsv, masterPathCsv, faceColours, outputPathMat1);
        catch err
            err
            close all;
        end
    end
end