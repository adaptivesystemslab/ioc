function [indExceedVeloLimit, banCount] = checkVeloQuality(useTime, ekfDQ, degOfInterest)
    % want to check the velocity for both the segment bounds and
    % non-segment bounds
    
    thresholdVelocity = 4; % absolute threshold for the normalized 
    offsets = 0; % to prevent edge cases from being missed
    
    % generate the boundaries to check
    % - all segment entries
    % - start of data to the first segment
    % - middle non-segment entries
    % - end of last segment to end of data
    [indToCheckStartSegment] = checkInd(useTime.startTimeInd-offsets, size(ekfDQ, 1));
    [indToCheckEndSegment]   = checkInd(useTime.endTimeInd+offsets, size(ekfDQ, 1));
    
    [rmsVeloSegment, peakVeloSegment]  = generateRMSVelo(ekfDQ, indToCheckStartSegment, indToCheckEndSegment, degOfInterest);
    [indExceedSegment, rmsVeloSegment, banCount] = kmeansBan(rmsVeloSegment, indToCheckStartSegment, indToCheckEndSegment);

    indToCheckStartArray = [1                               useTime.endTimeInd(1:end-1)-offsets    useTime.endTimeInd(end)-offsets];
    indToCheckEndArray   = [useTime.startTimeInd(1)+offsets useTime.startTimeInd(2:end)+offsets    size(ekfDQ, 1)];
    [indToCheckStartNonSegment] = checkInd(indToCheckStartArray, size(ekfDQ, 1));
    [indToCheckEndNonSegment]   = checkInd(indToCheckEndArray, size(ekfDQ, 1));
    
    [rmsVeloNonSegment, peakVeloNonSegment] = generateRMSVelo(ekfDQ, indToCheckStartNonSegment, indToCheckEndNonSegment, degOfInterest);
%     [indExceedNonSegment, srmsNonSegment] = kmeansBan(srmsNonSegment, indToCheckStartNonSegment, indToCheckEndNonSegment); 
    indExceedNonSegment = [];% don't check the kmeans for non-segments

    % combine the two
    indExceedVeloLimit = [indExceedSegment indExceedNonSegment];
    srms = [rmsVeloSegment rmsVeloNonSegment];
    peakVelo = [peakVeloSegment peakVeloNonSegment];
    indToCheckStart = [indToCheckStartSegment indToCheckStartNonSegment];
    indToCheckEnd =   [indToCheckEndSegment indToCheckEndNonSegment];
    
    % but just in case, we'll also hard limit threshold it
    % add a catch for large spikes. Anything higher than 5
    % rad/s is pretty unbelievable...    
%     indToPurge = find(srms > thresholdVelocity);
%     for ind_ban = 1:length(indToPurge)
%         currIndToPurge = indToPurge(ind_ban);
%         indExceedVeloLimit = [indExceedVeloLimit indToCheckStart(currIndToPurge):indToCheckEnd(currIndToPurge)];
%         srms(currIndToPurge) = 0;
%         banCount = banCount + 1;
%     end   

    indToPurge = find(peakVelo > thresholdVelocity);
    for ind_ban = 1:length(indToPurge)
        currIndToPurge = indToPurge(ind_ban);
        indExceedVeloLimit = [indExceedVeloLimit indToCheckStart(currIndToPurge):indToCheckEnd(currIndToPurge)];
        srms(currIndToPurge) = 0;
        banCount = banCount + 1;
    end
    
    indExceedVeloLimit = unique(sort(indExceedVeloLimit));
    
    if max(max(ekfDQ)) > 70 &&  max(max(ekfDQ)) < 80
        blah = 1;
    end
end

function arrayToCheck = checkInd(arrayToCheck, maxLength)
    tooSmall = find(arrayToCheck < 1);
    arrayToCheck(tooSmall) = 1;
    
    tooLarge = find(arrayToCheck > maxLength);
    arrayToCheck(tooLarge) = maxLength;
end

function [rmsVelo, peakVelo]  = generateRMSVelo(ekfDQ, indToCheckStart, indToCheckEnd, degOfInterest)
    rmsVelo = zeros(size(indToCheckStart));
    peakVelo = zeros(size(indToCheckStart));
    
    for ind_segments = 1:length(indToCheckStart)
        if indToCheckEnd(ind_segments) < indToCheckStart(ind_segments)
            fprintf('generateRMSVelo: Segmentation indexing issue (dt %d=%ds)\n', ind_segments, indToCheckEnd(ind_segments) - indToCheckStart(ind_segments));
            continue
        end
        
        % calculate the mean square velocity of each segment
        currVelo = ekfDQ(indToCheckStart(ind_segments):indToCheckEnd(ind_segments), degOfInterest);
        rmsVelo(ind_segments) = sum(rms(currVelo)) / length(degOfInterest);
        peakVelo(ind_segments) = max(max(abs(currVelo)));
    end
end

function [indExceedVeloLimit, srms, banCount] = kmeansBan(srms, indToCheckStart, indToCheckEnd)
    % if there is 1 or 2 in one cluster, then it's probably
    % problematic...
    indExceedVeloLimit = [];
    banCount = 0;
    
    srms = srms';
    
    if length(srms) < 2
        % do nothing
    else
        [idx, c] = kmeans(srms, 2);
        [c, maxind] = max(c);
        peakInd = idx == maxind;
        
        if sum(peakInd) < 2
            indToPurge = find(peakInd);
            for ind_ban = 1:length(indToPurge)
                currIndToPurge = indToPurge(ind_ban);
                indExceedVeloLimit = [indExceedVeloLimit indToCheckStart(currIndToPurge):indToCheckEnd(currIndToPurge)];
                srms(currIndToPurge) = 0;
                banCount = banCount + 1;
            end
        end
    end
    
    srms = srms';
end