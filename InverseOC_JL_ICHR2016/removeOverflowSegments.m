function segmentInfo = removeOverflowSegments(ekfTimeRaw, segmentInfo)
% modify the segmentstokeep to remove overflows. if both the start and end
% of a segment is not within te ekftimetokeep, then remove that segment
    segmentsToKeep = 1:length(segmentInfo.timeStart);
   
    for i = 1:length(segmentsToKeep)
        
        test(1) = sum(segmentInfo.timeStart(segmentsToKeep(i)) >= ekfTimeRaw);
        test(2) = sum(segmentInfo.timeStart(segmentsToKeep(i)) <= ekfTimeRaw);
        test(3) = sum(segmentInfo.timeEnd(segmentsToKeep(i)) >= ekfTimeRaw);
        test(4) = sum(segmentInfo.timeEnd(segmentsToKeep(i)) <= ekfTimeRaw);
        
        if ~sum(test == 0) % keep it
            
        else % remove it
            segmentsToKeep(i) = 0;
        end
    end
    
    segmentsToKeep = segmentsToKeep(segmentsToKeep > 0);
    
    % if nothing passed, then just leave it
    if ~isempty(segmentsToKeep)
        segmentInfo.use = segmentInfo.use(segmentsToKeep);
        segmentInfo.segmentCount = segmentInfo.segmentCount(segmentsToKeep);
        segmentInfo.timeStart = segmentInfo.timeStart(segmentsToKeep);
        segmentInfo.timeEnd = segmentInfo.timeEnd(segmentsToKeep);
    end
end