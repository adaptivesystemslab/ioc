function [withinSegmentCheck] = checkWithinSegment(feature_full, t_windowCount, segmentInfo)
% checkStandingVSMovement(feature_full, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2), ...
% segmentInfo, ratioRMSE{ind_windowCount}, ratioAll_end, ratioAll_all)

withinSegmentCheck = zeros(size(feature_full.t));
% check to see if the time point is within the segment
for ind_time = 1:length(t_windowCount)
    startCheck = segmentInfo.timeStart <= t_windowCount(ind_time);
    endCheck =                            t_windowCount(ind_time) <= segmentInfo.timeEnd;
    unionCheck = and(startCheck, endCheck);
    checkPass = find(unionCheck);
    
    if isempty(checkPass)
        % this data point is not within a segment
        checkPass = 0;
    else
        checkPass = 1;
    end
 
    withinSegmentCheck(ind_time) = checkPass;
end