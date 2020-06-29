function loadAnalysis()
    basePath = 'D:\aslab_gitlab\expressive-ioc\Data\IOC\';
    masterPathCsv = 'D:\aslab_gitlab\expressive-ioc\Data\IOC\summary.csv';
    
    % iterate through all the folders in bathpath to produce figures and
    % table summary of the results
    
    dirPath = dir(basePath);
    
    for i = 1:length(dirPath)
        currPath = fullfile(basePath, dirPath(i).name);
        dirCurrPath = dir(fullfile(currPath, '*.mat'));
        
        for j = 1:length(dirCurrPath)
            matPath = fullfile(currPath, dirCurrPath(j).name);
            matDataTemp = load(matPath);
            
            if isfield(matDataTemp, 'saveVar')
                matData = matDataTemp.saveVar;
                outputPathFig1 = sprintf("%s_%s_%s_fig_cum", currPath, ['Subject' num2str(matData.subjectId)], '');
                outputPathFig2 = sprintf("%s_%s_%s_fig_ind", currPath, ['Subject' num2str(matData.subjectId)], '');
                outputPathCsv = sprintf("%s_%s_%s_csv.csv", currPath, ['Subject' num2str(matData.subjectId)], '');
            elseif isfield(matDataTemp, 'outputVar')
                matData = matDataTemp.outputVar;
                outputPathFig1 = sprintf("%s_%s_%s_fig_cum", currPath, matData.trialInfo.name, matData.trialInfo.model);
                outputPathFig2 = sprintf("%s_%s_%s_fig_ind", currPath, matData.trialInfo.name, matData.trialInfo.model);
                outputPathCsv = sprintf("%s_%s_%s_csv.csv", currPath,  matData.trialInfo.name, matData.trialInfo.model);
            else
                continue;
            end
                        
%             try
                plotting2(matData, outputPathFig1, outputPathCsv, masterPathCsv);
                plotting3(matData, outputPathFig2, outputPathCsv, masterPathCsv);
               
%             catch err 
%                 err
%                 close all;
%             end
         
        end
    end
end

function h = plotting2(matData, outputPathFig, outputPathCsv, masterPathCsv)
    % load and process data
    dataIndsRan = matData.frameInds;
    t = matData.t(dataIndsRan);
    q = matData.q(dataIndsRan, :);
    progressVar = matData.progress;
    
    % the weight at each timestep is the sum of every window that overlaps with iterate
    for i = 1:length(t)
        weightAtI = [];
        for j = 1:length(progressVar)
            if sum(find(progressVar(i).winInds == j))
                weightAtI = [weightAtI; progressVar(j).weights];
            end
        end
        
        weights(i, :) = mean(weightAtI, 1);
        winCount(i, :) = size(weightAtI, 1);
    end
    
    weightLabels = matData.featureLabels;
    
    % make fig
    h = figure; 
    ax(1) = subplot(311);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Averaged weight windows');
    
    ax(2) = subplot(312);
    area(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(313);
    area(t, winCount);
    ylabel('WinCountUsed');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
end

function h = plotting3(matData, outputPathFig, outputPathCsv, masterPathCsv)
    % load and process data
    dataIndsRan = matData.frameInds;
    t = matData.t(dataIndsRan);
    dt = matData.dt;
    q = matData.q(dataIndsRan, :);
    progressVar = matData.progress;
    minLenThres = matData.minLenThres;
    maxLenThres = matData.maxLenThres;
    minRankThres = matData.minRankThres;
    
    for i = 1:length(t)
        weights(i, :) = progressVar(i).weights;
        rank(i, :) = progressVar(i).rankTraj(end);
        winLen(i, :) = length(progressVar(i).winInds);
    end
    
    weightLabels = matData.featureLabels;
    
    % make fig
    h = figure; 
    ax(1) = subplot(411);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Individual weight windows');
    
    ax(2) = subplot(412);
    area(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(413);
    area(t, rank);
    hold on;
    plot([t(1) t(end)], [minRankThres minRankThres], 'k', 'LineWidth', 2);
    ylim([0 minRankThres*1.2]);
    ylabel('Final Rank');
    
    ax(4) = subplot(414);
    area(t, winLen*dt);
    hold on;
    plot([t(1) t(end)], [minLenThres minLenThres]*dt, 'k', 'LineWidth', 2);
    plot([t(1) t(end)], [maxLenThres maxLenThres]*dt, 'k', 'LineWidth', 2);
    ylim([minLenThres/1.2*dt maxLenThres*1.2*dt]);
    ylabel('Win Len [s]');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
end

function h = plotting(matData, outputPathFig, outputPathCsv, masterPathCsv)
    % load and process data
    t = matData.t;
    q = matData.q;
    weights = matData.weightTraj;
    rankRatio = matData.rankTraj;
    segArray = matData.segmentArray(1:end-1, :);
    weightLabels = matData.featureLabels;
    gamma = matData.runParam.gamma;
    
    if ~isempty(segArray)
        for i = 1:size(segArray, 1)
            segmentInfo(i).timeStart = t(segArray(i, 1));
            segmentInfo(i).timeEnd = t(segArray(i, 2));
            segmentInfo(i).timeWidth = segmentInfo(i).timeEnd - segmentInfo(i).timeStart;
            
            tTimeSeries(i) = segmentInfo(i).timeStart;
            weightTimeSeries(i, :) = weights(segArray(i, 2), :);
        end
        
        ind = size(segArray, 1);
        tTimeSeries(ind + 1) = segmentInfo(ind).timeEnd;
        weightTimeSeries(ind + 1, :) = weights(segArray(ind, 2), :);
        
        timeWidth = [segmentInfo.timeWidth];
    else
%         tTimeSeries(1) = 0;
%         tTimeSeries(2) = t(end);
%         weightTimeSeries(1, :) = zeros(size(weights, 2), 1);
%         weightTimeSeries(2, :) = zeros(size(weights, 2), 1);
    end
   
    % make fig
    h = figure; 
    ax(1) = subplot(411);
    plot(t, q);
    if ~isempty(segArray)
        plotBoxes(segmentInfo, 'k');
    end
    ylabel('Joint Angle [rad]');
    
    ax(2) = subplot(412);
    if ~isempty(segArray)
        area(tTimeSeries, weightTimeSeries);
        lgd = legend(weightLabels,'AutoUpdate','off');
        lgd.NumColumns = 1;
        plotBoxes(segmentInfo, 'k', 0, -0.05, 1.05);
    end
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(413);
    plot(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    if ~isempty(segArray)
        plotBoxes(segmentInfo, 'k', 0, -0.05, 1.05);
    end
    ylabel('Weights [0:1]');
    
    ax(4) = subplot(414);
    plot(t, rankRatio); 
    if ~isempty(segArray)
        ylim([-1 gamma]);
        plotBoxes(segmentInfo, 'k', 0, -10, gamma+10);
    end
    ylabel('Rank Ratio');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
end
