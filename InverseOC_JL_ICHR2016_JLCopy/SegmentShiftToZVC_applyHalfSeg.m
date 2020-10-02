function useTime = SegmentShiftToZVC_applyHalfSeg(manualSegmentData, jointAngleData, settings)
    % given some existing time series segmentation data and joint angles,
    % provide the manual segments, slides it to the closest previous ZVC if
    % within some threshold of it. also looks for ZVCs within the window
    % and declare half segments is required
    
    % unit of time in seconds
    
    % variables
    offsetVal = settings.offsetVal; % initial offset applied to the manual segments (-0.2 for healthy1)
    toleranceGap = settings.toleranceGap; % how far back should the ZVC be before we ignore it (0.5 for helathy1)
%     sigOfConsideration = [1 4]; % the DOFs we care about
    plotData = 0; % make a plot of the data after we're done so we can chcek the results
    segmentProfileMod = settings.segmentProfileMod; % split the primitives into half primitives
    
    offsetValHalfSegment = 0.0; % shifting for the half segment section, to narrow the "search" for the mid ZVC
    
    % figure out the sigDof
%     veloTemp = diff(jointAngleData.jointAngles(sigOfConsideration, :)');
    
    [sigVelo, sigDofInd] = findSigDofVelo(jointAngleData.jointVelo');
    jointAngleData.sigJoint = jointAngleData.jointAngles(sigDofInd, :);
    jointAngleData.sigVelo = sigVelo;
    
    % pull out all the instances where ZVCs occur
    [crossingStruct, crossingDof] = zvcCheck(jointAngleData.time, sigVelo);
    crossingStruct = crossingStruct{1};
    
    % verify that all segment points are within the boundaries of the
    % available data. discard all pairs that are not
%     zeroOutKeep = [];
%     for i = 1:length(manualSegmentData.timeStart)
%         if manualSegmentData.timeStart(i) < jointAngleData.time(1) || ...
%                 manualSegmentData.timeEnd(i) < jointAngleData.time(1) || ...
%                 manualSegmentData.timeStart(i) > jointAngleData.time(end) || ...
%                 manualSegmentData.timeEnd(i) > jointAngleData.time(end)
%             % get rid of this
%         else
%             zeroOutKeep = [zeroOutKeep i];
%         end
%     end
    
%     manualSegmentData.timeStart = manualSegmentData.timeStart(zeroOutKeep);
%     manualSegmentData.timeEnd = manualSegmentData.timeEnd(zeroOutKeep);
%     manualSegmentData.segmentName = manualSegmentData.segmentName(zeroOutKeep);
%     manualSegmentData.segmentIncludeInd = manualSegmentData.segmentIncludeInd(zeroOutKeep);
%     manualSegmentData.use = manualSegmentData.use(zeroOutKeep);
%     manualSegmentData.segmentCount = 1:length(zeroOutKeep);
    
       
    % pull out the entries
    for i = 1:length(manualSegmentData.timeStart)
        currTimeStart(i) = manualSegmentData.timeStart(i);
        currTimeEnd(i) = manualSegmentData.timeEnd(i);
%        currSegmentName{i} = manualSegmentData.segmentName{i};
        
        % shift starting points, adding an offset
        [useTimeStartVal(i), useTimeStartInd(i)] = ...
            findZVCMatch(currTimeStart(i) + offsetVal, crossingStruct, jointAngleData, toleranceGap);
        
        % shift ending points
        [useTimeEndVal(i), useTimeEndInd(i)] = ...
            findZVCMatch(currTimeEnd(i) + offsetVal, crossingStruct, jointAngleData, toleranceGap);
    end
    
    switch segmentProfileMod
        case 'merge'
            % create merged segment points
            existingLength = length(currTimeStart);

            newStartPtInd = 1:2:existingLength;
            newEndPtInd = 2:2:existingLength;
            useTimeStartVal = useTimeStartVal(newStartPtInd);
            useTimeStartInd = useTimeStartInd(newStartPtInd);
            useTimeEndVal   = useTimeEndVal(newEndPtInd);
            useTimeEndInd   = useTimeEndInd(newEndPtInd);
            
        case 'segmentId'
           % merge based on segment id
           uniqueId = unique(manualSegmentData.segmentId);
           uniqueId = uniqueId(uniqueId > 0); % remove the no-use ones
           
           for i_uniqueId = 1:length(uniqueId)
                allInd = find(manualSegmentData.segmentId == uniqueId(i_uniqueId));
                newStartPtInd(i_uniqueId) = allInd(1);
                newEndPtInd(i_uniqueId) = allInd(end);
           end
           
           useTimeStartVal = useTimeStartVal(newStartPtInd);
           useTimeStartInd = useTimeStartInd(newStartPtInd);
           useTimeEndVal   = useTimeEndVal(newEndPtInd);
           useTimeEndInd   = useTimeEndInd(newEndPtInd);
           
        case 'normal'
             
        case 'half'
            % create half segment points
            existingLength = length(currTimeStart);
            
            for i = 1:length(useTimeStartVal)
                currTimeStartHalfSeg(i) = useTimeStartVal(i) + offsetValHalfSegment;
                currTimeEndHalfSeg(i) = useTimeEndVal(i) - offsetValHalfSegment;
                
                % find mid-ZVC point
                [newEndPt(i), newStartPt(i)] = ...
                    findZVCMid(currTimeStartHalfSeg(i), currTimeEndHalfSeg(i), crossingStruct, jointAngleData, toleranceGap);
            end
            
            % now add it to a start and an end
            useTimeStartVal = sort([useTimeStartVal newStartPt]);
            useTimeEndVal = sort([useTimeEndVal newEndPt]);
            
            [~, useTimeStartInd] = findClosestValue(useTimeStartVal, jointAngleData.time);
            [~, useTimeEndInd] = findClosestValue(useTimeEndVal, jointAngleData.time);
    end
    
    % remove repetitions in the values, as well as segment lengths that
    % equals 0
    zeroLengthArray = useTimeStartInd == useTimeEndInd;
    if ~isempty(find(zeroLengthArray, 1))
        arrayToKeep = not(zeroLengthArray);
        useTimeStartVal = useTimeStartVal(arrayToKeep);
        useTimeStartInd = useTimeStartInd(arrayToKeep);
        useTimeEndVal = useTimeEndVal(arrayToKeep);
        useTimeEndInd = useTimeEndInd(arrayToKeep);
    end
    
    useTime.startTimeVal = useTimeStartVal;
    useTime.endTimeVal = useTimeEndVal;
    useTime.startTimeInd = useTimeStartInd;
    useTime.endTimeInd = useTimeEndInd;    
    useTime.segmentCount = 1:length(useTimeStartVal);
    useTime.segmentInclude = ones(size(useTime.segmentCount)); % manualSegmentData.segmentIncludeInd; TODO FIX THIS
    
    if plotData
        h = figure;
        hold on
        ax(1) = subplot(211);
        plot(jointAngleData.time, jointAngleData.jointAngles(:, :));
%         plotXO(h, currTimeStart, currTimeEnd, useTime.timeStartVal, useTime.timeEndVal);
        plotBoxes(h, currTimeStart, currTimeEnd, 'r', 0.1);
        plotBoxes(h, useTime.startTimeVal, useTime.endTimeVal, 'b', -0.2);
        
        ax(2) = subplot(212);
        plot(jointAngleData.time, jointAngleData.jointVelo);
%         plotXO(h, currTimeStart, currTimeEnd, useTime.timeStartVal, useTime.timeEndVal);
        plotBoxes(h, currTimeStart, currTimeEnd, 'r', 0.1);
        plotBoxes(h, useTime.startTimeVal, useTime.endTimeVal, 'b', -0.2);
  
        linkaxes(ax, 'x');
        title('r = orig , g = corr')
        
        close(h);
    end
end

function plotXO(h, currTimeStart, currTimeEnd, useTimeStartVal, useTimeEndVal)
    figure(h);
    
    hold on
    scatter(currTimeStart, zeros(size(currTimeStart)), 'rx'); % orig start
    scatter(currTimeEnd, zeros(size(currTimeEnd)), 'ro'); % orig end

    scatter(useTimeStartVal, zeros(size(useTimeStartVal)), 'gx'); % new start
    scatter(useTimeEndVal, zeros(size(useTimeStartVal)), 'go'); % new end

    title('rx = orig start, go = corr end')
    hold off
end

function [useTimeStartVal, useTimeStartInd] = findZVCMatch(currTime, crossingStruct, jointAngleData, toleranceGap)
    %         [closeVal, closeInd] = findClosestValue(valueToFind, targetArray, favour)
    [closeVal, closeInd] = findClosestValue(currTime, crossingStruct.Time);
    offsetAmount = currTime - closeVal;

    if abs(offsetAmount) > toleranceGap
        % too far away from tolerance gap, don't shift the points 
        [closeVal, closeInd] = findClosestValue(currTime, jointAngleData.time, 'below');

        useTimeStartInd =  closeInd;
        useTimeStartVal =  closeVal;
    else
        % return the ZVC point
        targetVal = crossingStruct.Time(closeInd);
        [closeVal, closeInd] = findClosestValue(targetVal, jointAngleData.time, 'below');
         
        useTimeStartInd =  closeInd;
        useTimeStartVal =  closeVal;
    end
end

function [closeVal1, closeVal2] = ...
            findZVCMid(currTimeStart, currTimeEnd, crossingStruct, jointAngleData, toleranceGap)
    
    midVal = 0.5*(currTimeStart + currTimeEnd);
    lowerBound = currTimeStart < crossingStruct.Time;
    upperBound = crossingStruct.Time < currTimeEnd;
    inBound = and(lowerBound, upperBound);
    minLengthThreshold = 5; % all segments should be a minimum length long
    rangePercentage = 0.9; % want to maintain this much joint range

    if sum(inBound)
        % if there is a value
        crossingStructIndexInBound = crossingStruct.Index(inBound);
        crossingStructTimeInBound = crossingStruct.Time(inBound);
        
        [~, startInd] = findClosestValue(currTimeStart, jointAngleData.time);
        [~, endInd] = findClosestValue(currTimeEnd, jointAngleData.time);
        
        [midVal, midInd] = findClosestValue(midVal, crossingStruct.Time(inBound)); % pick the one closest to the middle
        [~, midInd] = findClosestValue(midVal, jointAngleData.time);
        
        % pull out the max values at the middle to the start and the end.
        % we want to find the smallest zvc boundaries that would maintain
        % the peak velos as found at the middle...this is assuming that the
        % middle is a good segment point
        
%         startIndCheck = startInd+1:closeInd;
        startMaxVelo = max(abs(jointAngleData.sigVelo(startInd:midInd))); % max velo between starting pt and the middle
        startJointRange = range(jointAngleData.sigJoint(startInd:midInd))*rangePercentage;
        for i = length(crossingStructIndexInBound):-1:1
            currentMaxVelo = max(abs(jointAngleData.sigVelo(startInd:crossingStructIndexInBound(i))));
            currentJointRange = range(jointAngleData.sigJoint(startInd:crossingStructIndexInBound(i)));
            
            % if the max velo is maintained, and we didn't exceed the
            % midind
            if startMaxVelo == currentMaxVelo && ...
                    currentJointRange >= startJointRange && ...
                    crossingStructIndexInBound(i) <= midInd %&& ...
         %           crossingStructIndexInBound(i) - startInd >= minLengthThreshold
                startIndToKeep = crossingStructIndexInBound(i);
            end
        end
        
        endMaxVelo = max(abs(jointAngleData.sigVelo(midInd:endInd))); % max velo between starting pt and the middle
        endJointRange = range(jointAngleData.sigJoint(midInd:endInd))*rangePercentage;
        for i = 1:length(crossingStructIndexInBound)
            currentMaxVelo = max(abs(jointAngleData.sigVelo(crossingStructIndexInBound(i):endInd)));
            currentJointRange = range(jointAngleData.sigJoint(crossingStructIndexInBound(i):endInd));
            
            % if the max velo is maintained, and we didn't exceed the
            % midind
            if endMaxVelo == currentMaxVelo && ...
                    currentJointRange >= endJointRange && ...
                    crossingStructIndexInBound(i) >= midInd %&& ...
                %    endInd - crossingStructIndexInBound(i) >= minLengthThreshold
                endIndToKeep = crossingStructIndexInBound(i);
            end
        end
        
        closeVal1 = jointAngleData.time(startIndToKeep);
        closeVal2 = jointAngleData.time(endIndToKeep);
        
%         [closeVal1, closeInd] = findClosestValue(startIndToKeep, jointAngleData.time);
%         [closeVal2, closeInd] = findClosestValue(endIndToKeep, jointAngleData.time);
    else
        % if no values are selected...
        [closeVal1, closeInd] = findClosestValue(currTimeStart, jointAngleData.time); % insert a value so it's not zero. worst comes to worst, it's a value already used
        [closeVal2, closeInd] = findClosestValue(currTimeEnd, jointAngleData.time); % insert a value so it's not zero. worst comes to worst, it's a value already used

    end
end