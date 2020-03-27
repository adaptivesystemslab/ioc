function IOCAnalysis()
    setPaths();
%     nowstr = datestr(now, 'yyyymmddHHMMSS');
    nowstr = '20200316_fatigueEdges';
    nowstr2 = '20200316_fatigueEdges';
      
    basePath = ['D:\results\fatigue_ioc03_weightsPattern\' nowstr '\'];
    searchString = 'mat_*_3DOF_3CF*.mat';
    outputPath = ['D:\results\fatigue_ioc04_weightsCluster\' nowstr2 '\'];
    checkMkdir(outputPath);
    
    
    currBasePathDir = dir([basePath searchString]);
    for j = 1:length(currBasePathDir)
        currFileName = currBasePathDir(j).name;
        currFullPath = fullfile(basePath, currFileName);
        
        if strcmpi(currFileName(end), '.')
            continue; % it's . or ..
        end

        % load each entry and prepare to plot
        load(currFullPath);
        
        subjectNum = str2num(trialInfo.runName(8:9));
        combinedStats = [stats_q stats_dq stats_tau stats_weights];
        cumStats{subjectNum} = combinedStats;
    end
    
    % now plot everything
    for i = 1:length(cumStats{1})
        h = figure('Position', [-1919 69 1920 964.8000]);
        hold on
        
        for j = 1:length(cumStats)
            currStats = cumStats{j}(i);
            plot(
        end
    end
end