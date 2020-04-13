function h = plotting_individual(matData, outputPathFig, outputPathCsv, masterPathCsv)
    % load and process data
    lenWeights = length(matData.featureLabels);
    
    progressVar = matData.progress;
    minLenThres = matData.minLenThres;
    maxLenThres = matData.maxLenThres;
    minRankThres = matData.minRankThres;
    
%     dataIndsRan = matData.frameInds;
    
%     if length(dataIndsRan) ~= length(progressVar)
% find start
    dataIndsRan = [];
    for i = 1:length(progressVar) 
        if ~isempty(progressVar(i).weights)
            dataIndsRan = [dataIndsRan; i];
        end
    end
%         dataIndsRan = 1:length(progressVar);
%     end
    
    t = matData.t(dataIndsRan);
    dt = matData.dt;
    q = matData.q(dataIndsRan, :);
    
    for i = 1:length(progressVar)
        if isempty(progressVar(i).weights)
            weights(i, :) = zeros(1, lenWeights);
            rank(i, :) = 0;
            winLen(i, :) = 0;
        else
            weights(i, :) = progressVar(i).weights;
            rank(i, :) = progressVar(i).rankTraj(end);
            winLen(i, :) = progressVar(i).winInds(end) - progressVar(i).winInds(1) + 1;
        end
    end
    
    weights = weights(dataIndsRan, :);
    rank = rank(dataIndsRan, :);
    winLen = winLen(dataIndsRan, :);
    weightLabels = matData.featureLabels;
    
    % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
%     ax(1) = subplot(311);
%     plot(t, q);
%     ylabel('Joint Angle [rad]');
%     title('Individual weight windows');
%     
    ax(1) = subplot(211);
    area(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
%     ax(3) = subplot(413);
%     area(t, rank);
%     hold on;
%     plot([t(1) t(end)], [minRankThres minRankThres], 'k', 'LineWidth', 2);
%     ylim([0 minRankThres*1.2]);
%     ylabel('Final Rank');
    
%     ax(4) = subplot(414);
%     area(t, winLen*dt);
%     hold on;
%     plot([t(1) t(end)], [minLenThres minLenThres]*dt, 'k', 'LineWidth', 2);
%     plot([t(1) t(end)], [maxLenThres maxLenThres]*dt, 'k', 'LineWidth', 2);
%     ylim([minLenThres/1.2*dt maxLenThres*1.2*dt]);
%     ylabel('Win Len [s]');
    
    
    data = cell2mat({progressVar.winInds}');
    minimum = data(:, 1);
    maximum = data(:, 2);
%     h = figure;
    
    ax(2) = subplot(212);
%     aH = axes;

    bH = bar(t, maximum);
    hold on;
    bH1 = bar(t, minimum);
    for ii = bH1
        ii.FaceColor = [1 1 1];
        ii.LineStyle = 'None';
        ii.EdgeColor = [1 1 1];
    end
    for ii = bH
        ii.LineStyle = 'None';
        ii.EdgeColor = [1 1 1];
    end
    ylim([min(minimum) max(maximum)]);
%     axes('Position',aH.Position,'XTick',[],'YTick',[],'Color','None');

   linkaxes(ax, 'x');
   xlim([t(1) t(end)]);
    xlabel('Time [s]');
    
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
end