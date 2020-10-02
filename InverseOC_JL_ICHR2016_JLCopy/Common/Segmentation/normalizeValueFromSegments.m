function [meanVal, maxInd1, maxInd2] = normalizeValueFromSegments(ekfQ, ekfTime, timeStart, timeEnd, segmentsToPullFrom)
% calculate a normalization factor, max(max(abs(ekfQ))), from the denoted
% segments. if segmentsToPullFrom is empty or not defined, pull from all
% the segments. Otherwise, it expects an array of the segments

if ~exist('segmentsToPullFrom', 'var') || isempty(segmentsToPullFrom)
    segmentsToPullFrom = 1:length(timeStart);
end

% pull out the max value from each segment
for i = 1:length(segmentsToPullFrom)
    i;
    [closeVal_start, closeInd_start] = findClosestValue(timeStart(i), ekfTime);
    [closeVal_end, closeInd_end] = findClosestValue(timeEnd(i), ekfTime);
    
    ekfIndToUseLocal = closeInd_start:closeInd_end;
    ekfToUseLocal = ekfQ(ekfIndToUseLocal, :);
    
    [normVal_array(i), maxDof_array(i), maxT_array(i)] = normalizeSegment(ekfToUseLocal);
    
%     ekfIndToUse = [ekfIndToUse; ekfIndToUseLocal];
%     ekfToUse = [ekfToUse; ekfToUseLocal];
end

[maxVal, maxInd] = max(normVal_array);
maxInd1 = maxDof_array(maxInd);
maxInd2 = maxT_array(maxInd);

meanVal = mean(normVal_array);

function [normVal, maxDof, maxTimestep] = normalizeSegment(ekfQ)
    absQ = abs(ekfQ);
    [maxVal1, maxTimesteps] = max(absQ, [], 1);
    [normVal, maxDof] = max(maxVal1);
    maxTimestep = maxTimesteps(maxDof);

