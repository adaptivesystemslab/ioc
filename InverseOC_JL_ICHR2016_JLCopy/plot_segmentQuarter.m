% produce quarter segments. assume the input data is in half segments
% already

segmentInfo_quarter{1} = []; % upward
segmentInfo_quarter{2} = []; % upward 
segmentInfo_quarter{3} = []; % resting
segmentInfo_quarter{4} = []; % downward
segmentInfo_quarter{5} = []; % downward

for ind_seg = 1:2:length(segmentInfo.timeStart)
    seg1Delta = segmentInfo.timeEnd(ind_seg) - segmentInfo.timeStart(ind_seg);
    seg2Delta = segmentInfo.timeEnd(ind_seg+1) - segmentInfo.timeStart(ind_seg+1);
     
    % determine the start/end ind of each segment
    startSeg(1) = segmentInfo.timeStart(ind_seg);
    endSeg(1)   = startSeg(1) + (seg1Delta/2);
    startSeg(2) = startSeg(1) + (seg1Delta/2);
    endSeg(2)   = segmentInfo.timeEnd(ind_seg);
    
    startSeg(3) = segmentInfo.timeEnd(ind_seg);
    endSeg(3)   = segmentInfo.timeStart(ind_seg+1);
    
    startSeg(4) = segmentInfo.timeStart(ind_seg+1);
    endSeg(4)   = startSeg(4) + (seg2Delta/2);
    startSeg(5) = startSeg(4) + (seg2Delta/2);
    endSeg(5)   = segmentInfo.timeEnd(ind_seg+1);
    
    % then copy out each segment's profile
    for ind_seg2 = 1:length(startSeg)
        [startVal, startInd] = findClosestValue(startSeg(ind_seg2), feature_full.t);
        [endVal, endInd] = findClosestValue(endSeg(ind_seg2), feature_full.t);
        
        newAdd = avgWeightArray_belowThres( startInd:endInd, :);
        segmentInfo_quarter{ind_seg2} = [segmentInfo_quarter{ind_seg2}; newAdd]; 
    end
end

% now calc mean/avg
for ind_seg2 = 1:length(startSeg)
    meanWeight_group(ind_seg2, :) = mean(segmentInfo_quarter{ind_seg2}, 1);
    stdWeight_group(ind_seg2, :) = std(segmentInfo_quarter{ind_seg2}, 1);
end

% pull out the top values
maxAvgWeight = max(avgWeightArray_belowThres);
[~, maxInd1] = max(maxAvgWeight);
maxAvgWeight(maxInd1) = 0;
[~, maxInd2] = max(maxAvgWeight);
maxAvgWeight(maxInd2) = 0;
[~, maxInd3] = max(maxAvgWeight);
maxAvgWeight(maxInd3) = 0;

top3ind = sort([maxInd1, maxInd2, maxInd3]);

h25 = figure;
barwitherr(stdWeight_group(:, top3ind), 1:5, meanWeight_group(:, top3ind));
title('Top 3 basis weights by 1/5 seg');
legend(cost_function_names(top3ind));
% feature_full.t, avgWeightArray_belowThres,