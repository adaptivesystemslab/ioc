function [segLabel, segName, segInclude, segBlacklist] = segmentPointWindowExpand(obj, jointDataStruct, segInfo, maxArrayLength, settings)
    % generates the segLabel (p0 or p1 marker), segName (the English name)
    % and segInclude (weither or not to include these points into the
    % training :
    %    '1' is 'use as segment', 
    %    '0' is 'use as not-segment', 
    %   '-1' if can't use)
    % for segBlacklist, this is blackout for the testing component
    
    allowUnknownMarkers = settings.allowUnknownMarkers;

    if isempty(obj.settings.segmentPointWindow)
        offset = 3; % sets the window of the 'segment' points around the actual 'segment' point
        offsetVeloRelative = 100; % no velocity threshold
        offsetVeloHard = 100; % no velocity threshold
    else
        offset = obj.settings.segmentPointWindow;
        offsetVeloRelative = settings.segmentPointWindowVeloThresholdRelative; 
        offsetVeloHard = settings.segmentPointWindowVeloThresholdHard; 
    end
    
    segmentPt = [];
    segmentName = [];
    segmentInclude = zeros(size(segmentPt));
      
    % generate the segment name brackets first, for both segment points we
    % want to keep and don't
    segName = cell(maxArrayLength, 1);
    for ind_nameArray = 1:length(segInfo.startTimeInd)
        for ind_localName = segInfo.startTimeInd(ind_nameArray):segInfo.endTimeInd(ind_nameArray)
            
            if isempty(segName{ind_localName})
                stackStr = {};
                stackStr{1} = segInfo.segmentName{ind_nameArray};
            else
                % need to account for overlapping motions...
                stackStr = segName{ind_localName};
                stackStr = [stackStr segInfo.segmentName{ind_nameArray}];
            end
            segName{ind_localName} = stackStr;
        end
    end


%     for ind_nameArray = 1:length(segInfo.startTimeInd)
%         newSegmentPts = segInfo.startTimeInd(ind_nameArray):segInfo.endTimeInd(ind_nameArray);
%         
%         % generate a new name array
%         newNameArray = cell(size(newSegmentPts));
%         for ind_newName = 1:length(newSegmentPts)
%             newNameArray(ind_newName) = segInfo.segmentName(ind_nameArray);
%         end
%         
%         segmentPt = [segmentPt newSegmentPts];
%         segmentName = [segmentName newNameArray];
%     end
    

    
    % generate the offset that we're adding to the segment points
    [segmentPtStart, segmentNameStart, segmentIncludeStart, segmentClusterStart, segmentClusterStart_StartEdge, segmentClusterStart_EndEdge] = ...
        generateInitOffset(jointDataStruct, segInfo.startTimeInd, segInfo.segmentName, segInfo.segmentInclude, offset, offsetVeloRelative, offsetVeloHard, maxArrayLength);
    [segmentPtEnd, segmentNameEnd, segmentIncludeEnd, segmentClusterEnd, segmentClusterEnd_StartEdge, segmentClusterEnd_EndEdge] = ...
        generateInitOffset(jointDataStruct, segInfo.endTimeInd, segInfo.segmentName, segInfo.segmentInclude, offset, offsetVeloRelative, offsetVeloHard,maxArrayLength);
    
    segmentPt = [segmentPt segmentPtStart segmentPtEnd];
    segmentName = [segmentName segmentNameStart segmentNameEnd];
    segmentInclude = [segmentInclude segmentIncludeStart segmentIncludeEnd];
    
    % modify the segInfo with the local minimized version 
    segInfo.startTimeInd = segmentClusterStart;
    segInfo.endTimeInd = segmentClusterEnd;
    segInfo.startTimeInd_StartEdge = segmentClusterStart_StartEdge;
    segInfo.startTimeInd_EndEdge = segmentClusterStart_EndEdge;
    segInfo.endTimeInd_StartEdge = segmentClusterEnd_StartEdge;
    segInfo.endTimeInd_EndEdge = segmentClusterEnd_EndEdge;
    
    % highlight all the segments we can use with a 0, to denote for the
    % training classifier later that we can use it. later in the function,
    % the 1 will be added, to denote the segment points
    segIncludeNotSeg = zeros(maxArrayLength, 1);
    for ind_nonSegmentHighlight = 1:length(segInfo.startTimeInd)
        if segInfo.segmentInclude(ind_nonSegmentHighlight) == 1
            startInd = segInfo.startTimeInd(ind_nonSegmentHighlight);
            endInd = segInfo.endTimeInd(ind_nonSegmentHighlight);
            segIncludeNotSeg(startInd:endInd) = 1;
        end
        % NOTE TODO - this does not currently account for spaces between
        % motions that could be used for training, in the event that
        % connectStartToEnd is not raised
    end
    
    % second pass, now check for the options, but don't overfill any of the
    % stuff that has already been filled
    if settings.connectStartToEnd        
        for ind_segPtTemp = 1:length(segInfo.startTimeInd)+1; % +1 to check the last entry properly
            % combine the end of one segment to the start of the next           
            
            if ind_segPtTemp == 1
                % first seg point, only has a beginning
                
                if ~settings.connectToDataStart
                    continue % no modification needed at the start
                end
                
                startPt = 1; % fill it to the start of the time series
%                 endPt = segInfo.startTimeInd(1)-offset-1;
                endPt = segInfo.startTimeInd_StartEdge(1)-1;
                
            elseif ind_segPtTemp == length(segInfo.startTimeInd)+1
                % last seg point, only has an end
           
                if ~settings.connectToDataStart
                    continue % no modification needed at the start
                end

%                 startPt = segInfo.endTimeInd(end)+offset+1;
                startPt = segInfo.endTimeInd_EndEdge(end)+1;
                endPt = maxArrayLength; % fill it to the start of the time series
                
                
            else
                % connect the end of the segment to the start of the
                % current one
%                 startPt = segInfo.endTimeInd(ind_segPtTemp-1)+offset+1;
%                 endPt = segInfo.startTimeInd(ind_segPtTemp)-offset-1;
                startPt = segInfo.endTimeInd_EndEdge(ind_segPtTemp-1)+1; % the 'offset' is now incorporated into the edge values
                endPt = segInfo.startTimeInd_StartEdge(ind_segPtTemp)-1;
            end
            
            if startPt < 1
                startPt = 1;
            end
            
            if endPt > maxArrayLength
                endPt = maxArrayLength;
            end
            
            newSegmentPts = startPt:1:endPt;
            
            % generate a new name array
            newNameArray = cell(size(newSegmentPts));
            for ind_newName = 1:length(newSegmentPts)
                if ind_segPtTemp == length(segInfo.startTimeInd)+1
                    % the last entry (end of last entry to the end of the 
                    % dataset) should take on the name of the template 
                    % ahead of it
                    newNameArray(ind_newName) = segInfo.segmentName(ind_segPtTemp - 1);
                else
                    newNameArray(ind_newName) = segInfo.segmentName(ind_segPtTemp);
                end
            end
            
            segmentPt = [segmentPt newSegmentPts];
            segmentName = [segmentName newNameArray];
            
            % now, generate the include array
            if ind_segPtTemp == length(segInfo.startTimeInd)+1
                % so settings.connectToDataStart is raised and 
                % is at the last entry (last entry to end of data)
                segmentInclude = [segmentInclude ones(size(newSegmentPts))];
            elseif segInfo.segmentInclude(ind_segPtTemp) == 1
                % generate the include array
                segmentInclude = [segmentInclude ones(size(newSegmentPts))];
            else
                % NOTE, both 0 and -1 entries would end up here
                if ind_segPtTemp <= length(segInfo.startTimeInd) % testing: was only "<" before. if something odd happens, change back to that
                    % skip entries that are not listed on the include list
                    segmentInclude = [segmentInclude zeros(size(newSegmentPts))];
                elseif ind_segPtTemp == length(segInfo.startTimeInd)+1
                    % last entry is not included, so ignore the last
                    % segment-to-end of data as well
                    segmentInclude = [segmentInclude zeros(size(newSegmentPts))];
                end
            end
            
            
            if ~(ind_segPtTemp == 1) && ~(ind_segPtTemp == length(segInfo.startTimeInd)+1)
                % so overlap regions between two exemplars
                % generate a new name array
                newNameArray = cell(size(newSegmentPts));
                for ind_newName = 1:length(newSegmentPts)
                    newNameArray(ind_newName) = segInfo.segmentName(ind_segPtTemp - 1);
                end
                
                segmentPt = [segmentPt newSegmentPts];
                segmentName = [segmentName newNameArray];
                
                if segInfo.segmentInclude(ind_segPtTemp) == 1
                    % generate the include array
                    segmentInclude = [segmentInclude ones(size(newSegmentPts))];
                else
                    % NOTE, both 0 and -1 entries would end up here
                    if ind_segPtTemp <= length(segInfo.startTimeInd) % testing: was only "<" before. if something odd happens, change back to that
                        % skip entries that are not listed on the include list
                        segmentInclude = [segmentInclude zeros(size(newSegmentPts))];
                    elseif ind_segPtTemp == length(segInfo.startTimeInd)+1
                        % last entry is not included, so ignore the last
                        % segment-to-end of data as well
                        segmentInclude = [segmentInclude zeros(size(newSegmentPts))];
                    end
                end
            end
        end
    end
    
    
    [segmentPtSort, segmentPtInd] = sort(segmentPt);
    [segmentPtUnique, IA, IC] = unique(segmentPtSort);
    
    % convert 2 points from both segment and non-segment edges to unknown
    if allowUnknownMarkers
        segmentPtGaps = find(diff(segmentPtUnique) > 1);
        segmentPtConvert = [];
        for ind_convert = 1:length(segmentPtGaps)
            currGapPtStart = segmentPtUnique(segmentPtGaps(ind_convert));
            currGapPtEnd = segmentPtUnique(segmentPtGaps(ind_convert)+1);
            currAdd1 = [currGapPtStart-1 currGapPtStart currGapPtStart+1 currGapPtStart+2];
            currAdd2 = [currGapPtEnd-2   currGapPtEnd-1 currGapPtEnd     currGapPtEnd+1];
            segmentPtConvert = [segmentPtConvert currAdd1 currAdd2];
        end
    else
        segmentPtConvert = [];
    end
    
    segLabel = ones(maxArrayLength, 1) * obj.label_notSegment;
    segLabel(segmentPtUnique) = obj.label_segment;
    segLabel(segmentPtConvert) = obj.label_unknown;
    
    % now update the segment names, and combine all the include data we've
    % generated so far
    segIncludeSeg = zeros(maxArrayLength, 1);
    for ind_uniqueArray = 1:length(segmentPtUnique)
        findInd = find(segmentPt == segmentPtUnique(ind_uniqueArray));
        
        stackStr = {};
        stackInclude = 0;
        for ind_found = 1:length(findInd)
            % stack the indices, if the length changes
            % if there is an existing motion in segName, the one holding
            % the primitive we're interested (ie here) trumps
            stackStr{ind_found} = segmentName{findInd(ind_found)};
            stackInclude = stackInclude + segmentInclude(findInd(ind_found)); % so segmentInclude only covers +1 points
        end
        
            % if the name that is currently in segName does not equal any of
            % the names in stackStr, then add that too
            overlapTestArray = [segName{segmentPtUnique(ind_uniqueArray)} stackStr];
            if length(unique(stackStr)) == length(stackStr) && length(stackStr) > 1
                % if stackStr contains 2 motions of the same type, chances
                % are the one in segName is also the same. But if stackStr
                % contains two of the same, then we want to keep them both.
                % Update segName to be stackStr (ie do nothing)
            elseif length(unique(overlapTestArray)) == length(overlapTestArray)
                % so we know that stackStr contains no repeated labels. If
                % things here is true, then stackStr and segName holds
                % different entries, in which case we should combine them
                % all together
                stackStr = overlapTestArray;
            else
                % else there is some overlap between stackStr and
                % segName, which we'll have to assume is a repeat
             
            end
        
        % update the segName and segInclude array with the additional
        % points, and include them properly if need be
        segName{segmentPtUnique(ind_uniqueArray)} = stackStr;
        segIncludeSeg(segmentPtUnique(ind_uniqueArray)) = stackInclude > 0; % if stackInclude is 1 or 2, it'll just set it to 1
    end
    
    % remove blacklisted items from training
    if settings.ignorePointsFromFirstSegment
        segIncludeSeg(1:segInfo.endTimeInd(1)) = 0;
        segIncludeNotSeg(1:segInfo.endTimeInd(1)) = 0;
    end
    
    % now combine the segIncludeSeg and segIncludeNotSeg array
    segInclude = ones(maxArrayLength, 1)*-1; % exclude from any tests
    segInclude(segIncludeNotSeg == 1) = 0; % non-segment points
    segInclude(segIncludeSeg == 1) = 1; % segment points
    
    % generating the don't include list and insure it is the last thing
    % written so it is on the top
    blacklistItems = find(segInfo.segmentInclude == -1);
    blacklistInd = [];
    for ind_blacklist = 1:length(blacklistItems)
        % generate the indices that would be blacked out
        blackoutInd = blacklistItems(ind_blacklist);
        startPt = segInfo.startTimeInd(blackoutInd);
        endPt = segInfo.endTimeInd(blackoutInd);
        
        if startPt < 1
            startPt = 1;
        end
        
        if endPt > maxArrayLength
            endPt = maxArrayLength;
        end

        blacklistInd = [blacklistInd startPt:endPt];
    end
    
    segBlacklist = ones(maxArrayLength, 1); % exclude from any tests
    segBlacklist(blacklistInd) = 0;
end

function [segmentPt, segmentName, segmentInclude, segmentInd, segmentIndStartEdge, segmentIndEndEdge] = ...
    generateInitOffset(jointDataStruct, segmentPtTemp, segmentNameTemp, includeArray, offset, offsetVeloRelative, offsetVeloHard, maxArrayLength)
    % generate the initial offset. this should be ran for both the starting
    % and the ending array
    segmentPt = []; % the inds of the segment to use
    segmentName = {};
    segmentInclude = []; % segment points to use in training or not
    
    segmentInd = []; % a single point from each cluster for later processing. derived from segmentPt
    segmentIndStartEdge = []; 
    segmentIndEndEdge = [];
    
    [sigVelo, sigDofInd] = findSigDofVelo(jointDataStruct.jointVelo);
    sigVelo = sigVelo';
    
    for ind_segPtTemp = 1:length(segmentPtTemp)
        
        currPt = segmentPtTemp(ind_segPtTemp);
        startPt = currPt-offset;
        endPt = currPt+offset;

        if startPt < 1
            startPt = 1;
        end

        if endPt > maxArrayLength
            endPt = maxArrayLength;
        end
        
        % include all points in the gap:
        newSegmentPts = startPt:1:endPt;
        
        % then calculate new points via velocity threshold
        veloToExamine = sigVelo(newSegmentPts, :); % use sig dof instead of just rms
        veloToExamineRMS = rms(veloToExamine, 2);
        veloLowInd = find(veloToExamineRMS < offsetVeloHard); % pull out all the values that are below some absolute threshold
        
        if length(veloLowInd) < 3
            % however, a hard threshold may not generate any segment 
            % points. if that is the case, we'll use a variable threshold 
            % instead of a hard one
            veloLowInd = find(veloToExamineRMS < min(veloToExamineRMS) * offsetVeloRelative); 
            [clusterStart, clusterEnd, clusterLength] = findDiscontinuity(veloLowInd, 3); % ignore the gap if it's only differing by 1
        else
            [clusterStart, clusterEnd, clusterLength] = findDiscontinuity(veloLowInd, 3); % ignore the gap if it's only differing by 1
            
            % remove short clusters
            longClusters = clusterLength > 1;
            clusterStart = clusterStart(longClusters);
            clusterEnd = clusterEnd(longClusters);
            clusterLength = clusterEnd - clusterStart;
        end
      
%         % now take the largest continious cluster
%         [largestClusterVal, largestClusterInd] = max(clusterLength);
        
        % or the cluster with the smallest rms
        clusterRMS = zeros(size(clusterStart));
        for ind = 1:length(clusterStart)
            clusterRMS(ind) = rms(veloToExamine(clusterStart(ind):clusterEnd(ind)));
        end
        [largestClusterVal, largestClusterInd] = min(clusterRMS);
        
        newSegmentPts = newSegmentPts(clusterStart(largestClusterInd):clusterEnd(largestClusterInd));
        
        % generate a new name array
        newNameArray = cell(size(newSegmentPts));
        for ind_newName = 1:length(newSegmentPts)
            newNameArray(ind_newName) = segmentNameTemp(ind_segPtTemp);
        end
        
        if isempty(newSegmentPts)
           x = 1; 
        end

        segmentPt = [segmentPt newSegmentPts];
        segmentName = [segmentName newNameArray];
        segmentIndStartEdge = [segmentIndStartEdge newSegmentPts(1)];
        segmentInd = [segmentInd floor(mean([newSegmentPts(1) newSegmentPts(end)]))];
        segmentIndEndEdge = [segmentIndEndEdge newSegmentPts(end)];
        
        if includeArray(ind_segPtTemp) == 1
            segmentInclude = [segmentInclude ones(size(newSegmentPts))];
        else
            segmentInclude = [segmentInclude zeros(size(newSegmentPts))];
        end
    end
end