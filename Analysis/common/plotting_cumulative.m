function [matSave] = plotting_cumulative(matData, outputPathFig_all, outputPathFig_rank, outputPathCsv, masterPathCsv, faceColours,outputPathMat1)
    outputPathFig_all_overall = [outputPathFig_all '_overall'];
    outputPathFig_all_variance = [outputPathFig_all '_variance'];
    
    % load and process data
    progressVar = matData.progress;
%     dataIndsRan = matData.frameInds;
    rankThres = matData.minRankThres;
    featureLabels = matData.featureLabels;
    
%     if length(dataIndsRan) ~= length(progressVar)
    dataIndsRan = [];
    for i = 1:length(progressVar) 
        if ~isempty(progressVar(i).weights)
            dataIndsRan = [dataIndsRan; i];
        end
    end
%         dataIndsRan = 1:length(progressVar);
%     end
    
    t = matData.t(:);
    q = matData.q(:, :);

    % the weight at each timestep is the sum of every window that overlaps with iterate
    [weights_all, weights_all_var, winCount_all] = cumWeights(t, progressVar, featureLabels);
%     [weights_rank, weights_rank_var, winCount_rank] = cumWeights_rankFilt(t, progressVar, featureLabels, rankThres);
    
    t = matData.t(dataIndsRan);
    q = matData.q(dataIndsRan, :);
    weights_all = weights_all(dataIndsRan, :);
    weights_all_var = weights_all_var(dataIndsRan, :);
    winCount_all = winCount_all(dataIndsRan, :);
%     weights_rank = weights_rank(dataIndsRan, :);
%     weights_rank_var = weights_rank_var(dataIndsRan, :);
%     winCount_rank = winCount_rank(dataIndsRan, :);
%     
    weightLabels = matData.featureLabels;

    matSave.t = t;
    matSave.q = q;
    matSave.weights = weights_all;
    matSave.var = weights_all_var;
    matSave.winCount = winCount_all;
    save(outputPathMat1, 'matSave');
    
    
    % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
    ax(1) = subplot(311);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Averaged weight windows (All pass)');
    
    ax(2) = subplot(312);
    area(t, weights_all);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(313);
    area(t, winCount_all);
    ylabel('WinCountUsed');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig_all_overall, 'fig');
    saveas(h, outputPathFig_all_overall, 'png');
    close(h);
    
    % make fig
    plot_weightVar(weights_all, weights_all_var, t,q,faceColours, outputPathFig_all_variance, weightLabels);
     
    %         % make fig
    %     h = figure('Position', [488 342 560*2 420*2]);
    %     ax(1) = subplot(311);
    %     plot(t, q);
    %     ylabel('Joint Angle [rad]');
    %     title('Averaged weight windows (Rank threshold)');
    %
    %     ax(2) = subplot(312);
    %     area(t, weights_rank);
    %     lgd = legend(weightLabels,'AutoUpdate','off');
    %     lgd.NumColumns = 1;
    %     ylabel('Weights [0:1]');
    %
    %     ax(3) = subplot(313);
    %     area(t, winCount_rank);
    %     ylabel('WinCountUsed');
    %
%     linkaxes(ax, 'x');
%     xlabel('Time [s]');
%     
%     saveas(h, outputPathFig_rank, 'fig');
%     saveas(h, outputPathFig_rank, 'png');
%     close(h);
end