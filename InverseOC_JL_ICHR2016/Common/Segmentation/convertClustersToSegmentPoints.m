function [manualSegStruct, algSegStruct, h] = convertClustersToSegmentPoints(comparisonStruct, t_error, n_exp, plotFig)

    if ~exist('plotFig', 'var')
        plotFig = 0;
    end

    if ~exist('n_exp', 'var')
        n_exp = 25;
    end

    if ~exist('t_error', 'var')
        t_error = 0.2;
    end
    
    mergeSmallClustersFlag = 1;
    mergeSmallWindowsFlag = 0;
    
    % rename variable for consistency
    initTime = comparisonStruct.time(1);
    q = comparisonStruct.q;
    data = comparisonStruct.data;
    time = comparisonStruct.time - initTime;
    labelGroundTruth = comparisonStruct.labelGroundTruth;
    segTimeStartGroundTruth = comparisonStruct.segTimeStartGroundTruth - initTime;
    segTimeEndGroundTruth = comparisonStruct.segTimeEndGroundTruth - initTime;
    
    labelAlgorithm = comparisonStruct.labelAlgorithm;
    
    % format things in the way manualSegmentationComparison expects it
    manualSegStruct.manualSegTime = sort([segTimeStartGroundTruth segTimeEndGroundTruth]);
    manLabel = {};
    for i = 1:length(manualSegStruct.manualSegTime)
        manLabel{i} = '-';
    end
    manualSegStruct.label = manLabel;
    manualSegStruct.id = ones(size(manualSegStruct.manualSegTime));
    
    % also, setting up this variable in case we lose all the clusters and
    % need to return a blank index
    algSegStruct.segmentTime = [];
    algSegStruct.segmentLabel = [];
    algSegStruct.segmentId = [];

    if mergeSmallClustersFlag
    % check for label array issues and erase small clusters
    % this is the long way of doing this, maybe come back and fix up later
    activeLabel = labelAlgorithm(1); % the current label
    activeLabelStart = 1; % and its index
    clusterCount = 1;
    for i = 2:length(labelAlgorithm)
        currLabel = labelAlgorithm(i);
        if currLabel ~= activeLabel 
            % label change
            if i - activeLabelStart(clusterCount) < n_exp
                % and it's too short, and not the first or last cluster
                if clusterCount == 1 || i == length(labelAlgorithm)
                    % otherwise just update the current values
                    clusterCount = clusterCount + 1;
                    activeLabel = labelAlgorithm(i);
                    activeLabelStart(clusterCount) = i;
                else
                    labelAlgorithm(activeLabelStart(clusterCount):i-1) = abs(activeLabel - 1); % this should flip the indices
                    clusterCount = clusterCount - 1;
                    activeLabel = labelAlgorithm(i); % update the active label
                end
            else
                % otherwise just update the current values
                clusterCount = clusterCount + 1;
                activeLabel = labelAlgorithm(i);
                activeLabelStart(clusterCount) = i;
                
            end
            
%             clf; plot(comparisonStruct.labelAlgorithm, 'x'); hold on; 
%             plot(labelAlgorithm - 0.02, 'o'); xlim([i-150 i+150]); ylim([-0.5 1.5]);
%             plot([i i], [0 1], 'rx');
        end
    end
    end
    
%     figure
    
    % cluster similar points together
    labelAlgorithmPad = padarray(labelAlgorithm, 1, 'both'); % padding the array so that there will always be a clean starting and ending cluster
    diffLabel = [0; diff(labelAlgorithmPad)];
    clusterStartSegment = find(diffLabel == 1) - 1; % start of a '1' cluster (-1 to remove the prepadding value)
    clusterEndSegment = find(diffLabel == -1) - 1; % end of a '1' cluster
    
    if isempty(clusterStartSegment)
        % no clusters survived the smoothing. return empty array
        if plotFig % show figure
            h = plotSegFigure;
        else
            h = [];
        end
        
        return
    end
    
    % for clusters that are really short, erase them. Assume hick-ups from
    % the algorithm side
%     figure; plot(diffLabel*0.1, 'x'); ylim([-1.5 1.5]);
%     hold on
%     plot(labelAlgorithm*0.1 - 0.05, ['rx']);
    
    % generate the segment clusters
% %     clusterStartSegment = find(diffLabel == 1) - 1; % start of a '1' cluster (-1 to remove the prepadding value)
% %     clusterEndSegment = find(diffLabel == -1) - 1; % end of a '1' cluster
% % 
% %     if clusterEndSegment(end) > length(labelAlgorithm)
% %         % the last cluster ends at the end of the data
% %         clusterEndSegment(end) = length(labelAlgorithm);
% %     end
% %     
% %     % generate the non-segment clusters
% %     clusterStartNonsegment = clusterEndSegment(1:end-1)+1;
% %     clusterEndNonsegment = clusterStartSegment(2:end)-1;
% % 
% %     % filter out short clusters
% % %     diffClusterSegment = diffLabelMinus - diffLabelPlus;
% % %     diffClusterNonSegment = diffLabelPlus(2:end) - diffLabelMinus(1:end-1);
% %     
% %     clusterStart = sort([clusterStartSegment; clusterStartNonsegment]);
% %     clusterEnd = sort([clusterEndSegment; clusterEndNonsegment]); % (2:end-1)
% %     for i = 2:length(clusterStart)-1
% %        currDiffClusterLength = clusterEnd(i) - clusterStart(i);
% %        if currDiffClusterLength < n_exp
% %            % too small. flip the clusters
% %            
% %        end
% %     end
    
    
    diffCluster = clusterEndSegment(2:end-1) - clusterStartSegment(2:end-1);
    clustersToUseInd = [1; find(diffCluster > t_error)+1; length(clusterStartSegment)]; % always include the edge cases
    clusterStartSegmentCrop = clusterStartSegment(clustersToUseInd); 
    clusterEndSegmentCrop = clusterEndSegment(clustersToUseInd);
    diffClusterCrop = clusterEndSegmentCrop - clusterStartSegmentCrop;
    
    % so now these clusters denote start and end of arrays
    % first cluster would only be a start
    for ind_cluster = 1:length(diffClusterCrop)
        currDiffCluster = diffClusterCrop(ind_cluster);
        clusterStart = clusterStartSegmentCrop(ind_cluster);
        clusterEnd = clusterEndSegmentCrop(ind_cluster);
        
        % return the middle value of the cluster (for first and last
        % segment)
        clusterMid = clusterStart + floor(currDiffCluster/2);
        
        % return the first 1/3 and the 2/3 of the cluster (for all other
        % segments)
        clusterThirdFirst = clusterStart + ceil(n_exp);
        clusterThirdSecond = clusterEnd - floor(n_exp);
        
        if clusterThirdFirst > clusterThirdSecond % this if should prevent overlapping windows
            clusterThirdFirst = clusterStart + floor(currDiffCluster/3);
            clusterThirdSecond = clusterEnd - floor(currDiffCluster/3);
        end
        
        % and now for error checking
        if clusterThirdSecond > length(time)
            clusterThirdSecond = time;
        end
        
        if ind_cluster == 1
            % first cluster, so start of an entry
            segTimeStartAlgorithm(ind_cluster) = time(clusterThirdSecond);
        elseif ind_cluster == length(diffClusterCrop);
            % last cluster, so end of an entry
            segTimeEndAlgorithm(ind_cluster-1) = time(clusterThirdFirst);
        else
            segTimeEndAlgorithm(ind_cluster-1) = time(clusterThirdFirst);
            segTimeStartAlgorithm(ind_cluster) = time(clusterThirdSecond);
        end
    end
    
    % debugging
    segTimeStartAlgorithmOrig = segTimeStartAlgorithm;
    segTimeEndAlgorithmOrig = segTimeEndAlgorithm;
    
    if mergeSmallWindowsFlag
    checkWindowLengthFlag = 1;
    currInd = 1;
    if length(segTimeEndAlgorithm) > 1
        while checkWindowLengthFlag
            % look at the current window and the next window
            currSegWindow = segTimeEndAlgorithm(currInd) - segTimeStartAlgorithm(currInd);
            nextSegWindow = segTimeEndAlgorithm(currInd+1) - segTimeStartAlgorithm(currInd+1);
            
            % if both windows are too small, merge them
            if currSegWindow < windowLengthThreshold && nextSegWindow < windowLengthThreshold
                segTimeStartAlgorithm = segTimeStartAlgorithm([1:currInd currInd+2:end]); % remove the start of the next entry
                segTimeEndAlgorithm = segTimeEndAlgorithm([1:currInd-1 currInd+1:end]); % remove the end of the curr entry
            else
                % these entries are fine, continue
                currInd = currInd + 1;
            end
            
            if currInd >= length(segTimeStartAlgorithm)
                % got to the end of the array, can't combine anymore
                checkWindowLengthFlag = 0;
            end
        end
    end
    end
    
    if plotFig % show figure
        h = plotSegFigure;
    else
        h = [];
    end
    
    
    % make sure that no window is really that short
    diffSegwindow = segTimeEndAlgorithm - segTimeStartAlgorithm;

    % reformat the variables so that it fits the assessment function    
    if exist('segTimeStartAlgorithm', 'var')
        algSegStruct.segmentTime = sort([segTimeStartAlgorithm segTimeEndAlgorithm]);
        algLabel = {};
        for i = 1:length(algSegStruct.segmentTime)
            algLabel{i} = '-';
        end
        algSegStruct.segmentLabel = algLabel;
        algSegStruct.segmentId = ones(size(algSegStruct.segmentTime));
    end
    
    function h = plotSegFigure
        % now plot all the data
        
        colourGroundTruth = 'b'; % blue is the ground truth
        colourAlgorithm = 'r'; % red is the algorithm
        
        h = figure;
        hold on
        title('red = algorithm, blue = ground truth, o = original, x = used');
        plot(time, q(:, 1:5));
        plot(time, labelGroundTruth * 0.5, [colourGroundTruth 'x']);
        plot(time, comparisonStruct.labelAlgorithm * 0.5 -0.03, [colourAlgorithm 'o']);
        plot(time, labelAlgorithm * 0.5 -0.05, [colourAlgorithm 'x']);
        plotBoxes(h, segTimeStartGroundTruth, segTimeEndGroundTruth, colourGroundTruth);
        if exist('segTimeStartAlgorithm', 'var')
            plotBoxes(h, segTimeStartAlgorithm, segTimeEndAlgorithm, colourAlgorithm, -0.05);
        end
        
        
    end
end