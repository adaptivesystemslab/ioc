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
    h = figure('Position', [488 342 560*2 420*2]); 
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