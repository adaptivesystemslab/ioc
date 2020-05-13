function [ekfTimeToKeep, dt, segmentInfo, ekfTimeIndicesToKeep] = setTimeVector(ekfTimeRaw, segmentTempCrop, segmentTempSeg, sampingRate, currSubjInfo, shiftZVCSettings)
    
    % flag to determine if we should remove segments that are not within
    % the bounds of the timetokeep for currSubjInfo.exerciseCropLength = full
    removeSegmentOverflowFlag = shiftZVCSettings.removeSegmentOverflowFlag;
    cropLength = shiftZVCSettings.cropAllowLength;
    oddevenCropLength = 2/50;
    
    if sampingRate ~= 0
        ekfTime = [ekfTimeRaw(2):(1/sampingRate):ekfTimeRaw(end-1)]'; % offset applied to index to prevent extralation
    else
        ekfTime = ekfTimeRaw;
    end
    
    dt = mean(diff(ekfTime));
    
    % figure out where the crop length should be
    segmentCount = length(segmentTempSeg.segmentCount);
    cropMid = floor(segmentCount/2);
    switch currSubjInfo.exerciseCropLength
        case 'full'
            cropStart = 1;
            cropEnd = 1;
            
            if ~isempty(segmentTempCrop)
                timeArrayToKeepStart = segmentTempCrop.timeStart(cropStart);
                timeArrayToKeepEnd = segmentTempCrop.timeEnd(cropEnd);
            else
                timeArrayToKeepStart = ekfTime(1);
                timeArrayToKeepEnd = ekfTime(end);
            end
            
            segmentsToKeep = 1:length(segmentTempSeg.timeStart);
            
        case 'autocrop'
            useArray = find(segmentTempSeg.use == 1);
            timeArrayToKeepStart = segmentTempSeg.timeStart(useArray(1)) - cropLength;
            timeArrayToKeepEnd = segmentTempSeg.timeEnd(useArray(end)) + cropLength;
            
            % if there is manual cropping, respect the bounds given there
            if ~isempty(segmentTempCrop)
                if timeArrayToKeepStart < segmentTempCrop.timeStart(1)
                    timeArrayToKeepStart = segmentTempCrop.timeStart(1);
                end
                
                if timeArrayToKeepEnd > segmentTempCrop.timeEnd(1)
                    timeArrayToKeepEnd = segmentTempCrop.timeEnd(1);
                end
            end
            
            segmentsToKeep = 1:length(segmentTempSeg.timeStart);
            
        case 'firsthalf'
            cropStart = 1;
            cropEnd = cropMid;
            
            timeArrayToKeepStart = segmentTempSeg.timeStart(cropStart) - cropLength;
            timeArrayToKeepEnd = segmentTempSeg.timeEnd(cropEnd) + cropLength;
            
            segmentsToKeep = cropStart:cropEnd;
            
        case 'secondhalf'
            cropStart = cropMid+1; % since it's always getting floored on the previous half, we'll ceiling it here
            cropEnd = segmentCount;
            
            timeArrayToKeepStart = segmentTempSeg.timeStart(cropStart) - cropLength;
            timeArrayToKeepEnd = segmentTempSeg.timeEnd(cropEnd) + cropLength;
            
            segmentsToKeep = cropStart:cropEnd;
            
        case 'odd'
            counter = 0;
            if strcmpi(shiftZVCSettings.segmentProfileMod, 'half') % will be split into two later
                segmentsToKeep = 1:2:segmentCount;
                for i = segmentsToKeep
                    counter = counter + 1;
                    timeArrayToKeepStart(counter) = segmentTempSeg.timeStart(i) - oddevenCropLength;
                    timeArrayToKeepEnd(counter) = segmentTempSeg.timeEnd(i) + oddevenCropLength;
                end
            else
                segmentsToKeep = 1:4:segmentCount; % already split into two
                for i = segmentsToKeep
                    counter = counter + 1;
                    timeArrayToKeepStart(counter) = segmentTempSeg.timeStart(i) - oddevenCropLength;
                    timeArrayToKeepEnd(counter) = segmentTempSeg.timeEnd(i+1) + oddevenCropLength;
                end
            end

            
        case 'even'
            counter = 0;
            if strcmpi(shiftZVCSettings.segmentProfileMod, 'half') % will be split into two later
                segmentsToKeep = 2:2:segmentCount;
                for i = segmentsToKeep
                    counter = counter + 1;
                    timeArrayToKeepStart(counter) = segmentTempSeg.timeStart(i) - oddevenCropLength;
                    timeArrayToKeepEnd(counter) = segmentTempSeg.timeEnd(i) + oddevenCropLength;
                end
            else
                segmentsToKeep = 3:4:segmentCount;
                for i = segmentsToKeep
                    counter = counter + 1;
                    timeArrayToKeepStart(counter) = segmentTempSeg.timeStart(i) - oddevenCropLength;
                    timeArrayToKeepEnd(counter) = segmentTempSeg.timeEnd(i+1) + oddevenCropLength;
                end
            end
            
    end
    
    ekfTimeIndicesToKeep = [];
    
    for i = 1:length(timeArrayToKeepStart)
        cropStart = timeArrayToKeepStart(i);
        cropEnd = timeArrayToKeepEnd(i);
        
        if cropStart < ekfTime(1)
            cropStart = ekfTime(1);
        end
        
        if cropEnd > ekfTime(end)
            cropEnd = ekfTime(end);
        end
        
        [~, cropStartInd] = findClosestValue(cropStart, ekfTime);
        [~, cropEndInd]   = findClosestValue(cropEnd, ekfTime);
        
        ekfTimeIndicesToKeep = [ekfTimeIndicesToKeep cropStartInd:cropEndInd];
    end
    
    ekfTimeToKeep = ekfTime(ekfTimeIndicesToKeep);
    
    if removeSegmentOverflowFlag
        segmentsToKeep = removeOverflowSegments(ekfTimeToKeep, segmentTempSeg, segmentsToKeep);
    end
    
    if isa(segmentTempSeg, 'segmentationDataHandle')
        segmentInfo = segmentTempSeg.retainSegments(segmentsToKeep);
    else
        segmentInfo.use = segmentTempSeg.use(segmentsToKeep);
        segmentInfo.segmentCount = segmentTempSeg.segmentCount(segmentsToKeep);
        segmentInfo.timeStart = segmentTempSeg.timeStart(segmentsToKeep);
        segmentInfo.timeEnd = segmentTempSeg.timeEnd(segmentsToKeep);
    end
end

function segmentsToKeep = removeOverflowSegments(ekfTimeRaw, segmentTempSeg, segmentsToKeep)
% modify the segmentstokeep to remove overflows. if both the start and end
% of a segment is not within te ekftimetokeep, then remove that segment
    for i = 1:length(segmentsToKeep)
        
        test(1) = sum(segmentTempSeg.timeStart(segmentsToKeep(i)) >= ekfTimeRaw);
        test(2) = sum(segmentTempSeg.timeStart(segmentsToKeep(i)) <= ekfTimeRaw);
        test(3) = sum(segmentTempSeg.timeEnd(segmentsToKeep(i)) >= ekfTimeRaw);
        test(4) = sum(segmentTempSeg.timeEnd(segmentsToKeep(i)) <= ekfTimeRaw);
        
        if ~sum(test == 0) % keep it
            
        else % remove it
            segmentsToKeep(i) = 0;
        end
    end
    
    segmentsToKeep = segmentsToKeep(segmentsToKeep > 0);
end