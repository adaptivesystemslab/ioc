function [motionEnd, motionWhole, ratioAll_end, ratioAll_all] = checkStandingVSMovement(...
    feature_full, indToUse_subwindow, segmentInfo, ratioRMSE, ratioAll_end, ratioAll_all)
% checkStandingVSMovement(feature_full, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2), ...
% segmentInfo, ratioRMSE{ind_windowCount}, ratioAll_end, ratioAll_all)

t_windowCount = feature_full.t(:, indToUse_subwindow);

if isempty(ratioAll_end)
    ratioAll_end{1} = {};
    ratioAll_end{2} = {};
end

if isempty(ratioAll_all)
    ratioAll_all{1} = {};
    ratioAll_all{2} = {};
    ratioAll_all{3} = {};
end

checkPassArray = [];
% check to see if the time point is within the segment
for ind_time = 1:length(t_windowCount)
    startCheck = segmentInfo.timeStart <= t_windowCount(ind_time);
    endCheck =                            t_windowCount(ind_time) <= segmentInfo.timeEnd;
    unionCheck = and(startCheck, endCheck);
    checkPass = find(unionCheck);
    
    if isempty(checkPass)
        % this data point is not within a segment
        segmentIn = 0;
        checkPass = 0;
    else
        segmentIn = checkPass;
        checkPass = 1;
    end
 
    segmentInArray(ind_time) = segmentIn(1);
    checkPassArray(ind_time) = checkPass;
end

% consider the end of the window first
endItem = segmentInArray(end);
switch checkPassArray(end)
    case 0
        motionEnd = ['Not in motion: segment ' num2str(endItem)];
        ratioAll_end{1} = [ratioAll_end{1}; ratioRMSE];
        
    case 1
        motionEnd = ['In motion: segment ' num2str(endItem)];
        ratioAll_end{2} = [ratioAll_end{2}; ratioRMSE];
end

mostCommonItem = mode(segmentInArray); 
switch sum(checkPassArray)
    case 0
        motionWhole = ['Not in motion: segment ' num2str(mostCommonItem)];
        ratioAll_all{1} = [ratioAll_all{1}; ratioRMSE];
        
    case length(checkPassArray)
        motionWhole = ['In motion: segment ' num2str(mostCommonItem)];
        ratioAll_all{2} = [ratioAll_all{2}; ratioRMSE];
        
    otherwise
        motionWhole = ['Partial motion: segment ' num2str(mostCommonItem)];
        ratioAll_all{3} = [ratioAll_all{3}; ratioRMSE]; % everything else
end