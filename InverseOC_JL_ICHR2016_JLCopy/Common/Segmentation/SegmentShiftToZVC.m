function useTime = SegmentShiftToZVC(manualSegmentData, jointAngleData, sigOfConsideration, settings)
    % given some existing time series segmentation data and joint angles,
    % provide the manual segments, slides it to the closest previous ZVC if
    % within some threshold of it
    
    % unit of time in seconds
    
    % variables
    offsetVal = settings.offsetVal; % initial offset applied to the manual segments (-0.2 for healthy1)
    toleranceGap = settings.toleranceGap; % how far back should the ZVC be before we ignore it (0.5 for helathy1)
%     sigOfConsideration = [1 4]; % the DOFs we care about
    plotData = 0; % make a plot of the data after we're done so we can chcek the results
    halfSeg = settings.halfSeg; % split the primitives into half primitives
    
    offsetValHalfSegment = 0.0; % shifting for the half segment section, to narrow the "search" for the mid ZVC
    
    % figure out the sigDof
%     veloTemp = diff(jointAngleData.jointAngles(sigOfConsideration, :)');
    
    [sigVelo, sigDofInd] = findSigDofVelo(jointAngleData.jointVelo');
    jointAngleData.sigVelo = sigVelo;
    
    % pull out all the instances where ZVCs occur
    [crossingStruct, crossingDof] = zvcCheck(jointAngleData.time, sigVelo);
    crossingStruct = crossingStruct{1};
    
    % verify that all segment points are within the boundaries of the
    % available data. discard all pairs that are not
    zeroOutKeep = [];
    for i = 1:length(manualSegmentData.timeStart)
        if manualSegmentData.timeStart(i) < jointAngleData.time(1) || ...
                manualSegmentData.timeEnd(i) < jointAngleData.time(1) || ...
                manualSegmentData.timeStart(i) > jointAngleData.time(end) || ...
                manualSegmentData.timeEnd(i) > jointAngleData.time(end)
%             get rid of this
        else
            zeroOutKeep = [zeroOutKeep i];
        end
    end
    
    manualSegmentData.timeStart = manualSegmentData.timeStart(zeroOutKeep);
    manualSegmentData.timeEnd = manualSegmentData.timeEnd(zeroOutKeep);
    manualSegmentData.segmentName = manualSegmentData.segmentName(zeroOutKeep);
    manualSegmentData.segmentIncludeInd = manualSegmentData.segmentIncludeInd(zeroOutKeep);
    manualSegmentData.use = manualSegmentData.use(zeroOutKeep);
    manualSegmentData.segmentCount = 1:length(zeroOutKeep);
    
       
    % pull out the entries
    for i = 1:length(manualSegmentData.timeStart)
        currTimeStart(i) = manualSegmentData.timeStart(i);
        currTimeEnd(i) = manualSegmentData.timeEnd(i);
        currSegmentName{i} = manualSegmentData.segmentName{i};
        
        % shift starting points, adding an offset
        [useTimeStartVal(i), useTimeStartInd(i)] = ...
            findZVCMatch(currTimeStart(i) + offsetVal, crossingStruct, jointAngleData, toleranceGap);
        
        % shift ending points
        [useTimeEndVal(i), useTimeEndInd(i)] = ...
            findZVCMatch(currTimeEnd(i) + offsetVal, crossingStruct, jointAngleData, toleranceGap);
    end
    
%     if halfSeg
%         % create half segment points
%         existingLength = length(currTimeStart);
%         
%         for i = 1:length(useTimeStartVal)
%             currTimeStartHalfSeg(i) = useTimeStartVal(i) + offsetValHalfSegment;
%             currTimeEndHalfSeg(i) = useTimeEndVal(i) - offsetValHalfSegment;
%             
%             % find mid-ZVC point
%             [newEndPt(i), newStartPt(i), newEndInd(i), newStartInd(i)] = ...
%                 findZVCMid(currTimeStartHalfSeg(i), currTimeEndHalfSeg(i), crossingStruct, jointAngleData, toleranceGap);
%             
%             currSegmentName{existingLength+i} = manualSegmentData.segmentName{i};
%             manualSegmentData.segmentIncludeInd(existingLength+i) = 1;
%         end
%         
%         % now add it to a start and an end
%         useTimeEndVal = sort([useTimeEndVal newEndPt]);
%         useTimeStartVal = sort([useTimeStartVal newStartPt]);
%         
%         useTimeEndInd = sort([useTimeEndInd newEndInd]);
%         useTimeStartInd = sort([useTimeStartInd newStartInd]);
%     end
    
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
    
    [startTimeVal,IA,IC] = unique(useTimeStartVal);
%     useTime.startTimeVal = unique(useTimeStartVal);
%     useTime.startTimeInd = unique(useTimeStartInd);
%     useTime.endTimeVal = unique(useTimeEndVal);
%     useTime.endTimeInd = unique(useTimeEndInd);
    useTime.startTimeVal = startTimeVal;
    useTime.startTimeInd = useTimeStartInd(IA);
    useTime.endTimeVal = useTimeEndVal(IA);
    useTime.endTimeInd = useTimeEndInd(IA);
    
    useTime.segmentName = currSegmentName(IA);
    useTime.segmentInclude = manualSegmentData.segmentIncludeInd;
    
    if plotData
        h = figure;
        hold on
%         subplot(211);
        plot(jointAngleData.time, jointAngleData.jointAngles(:, :));
%         plot(jointAngleData.time, sigVelo);
        
        plotXO(h, currTimeStart, currTimeEnd, useTime.startTimeVal, useTime.endTimeVal);
%         plotBoxes(h, currTimeStart, currTimeEnd, useTime.startTimeVal, useTime.endTimeVal);
        
        plotBoxes(h, currTimeStart, currTimeEnd, 'r', -.1);
        plotBoxes(h, useTime.startTimeVal, useTime.endTimeVal, 'b', -0.1);
        
        subplot(212);
        plot(jointAngleData.time, velo);
         plotXO(h, currTimeStart, currTimeEnd, useTime.startTimeVal, useTime.endTimeVal);
   plotBoxes(h, currTimeStart, currTimeEnd, 'r', 0);
        plotBoxes(h, useTime.startTimeVal, useTime.endTimeVal, 'b', -0.1);
  
         
         
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

function [closeVal1, closeVal2, closeInd1, closeInd2] = ...
            findZVCMid(currTimeStart, currTimeEnd, crossingStruct, jointAngleData, toleranceGap)
    
        narrowFlag = 0;
        
    midVal = 0.5*(currTimeStart + currTimeEnd);
    lowerBound = currTimeStart < crossingStruct.Time;
    upperBound = crossingStruct.Time < currTimeEnd;
    inBound = and(lowerBound, upperBound);

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
        for i = length(crossingStructIndexInBound):-1:1
            currentMaxVelo = max(abs(jointAngleData.sigVelo(startInd:crossingStructIndexInBound(i))));
            
            % if the max velo is maintained, and we didn't exceed the
            % midind
            if narrowFlag
                if startMaxVelo == currentMaxVelo && crossingStructIndexInBound(i) <= midInd ...
                        && jointAngleData.time(midInd) - jointAngleData.time(startInd) > 0.75
                    startIndToKeep = crossingStructIndexInBound(i);
                end
            else
                startIndToKeep = midInd;
            end

        end
        
        endMaxVelo = max(abs(jointAngleData.sigVelo(midInd:endInd))); % max velo between starting pt and the middle
        for i = 1:length(crossingStructIndexInBound)
            currentMaxVelo = max(abs(jointAngleData.sigVelo(crossingStructIndexInBound(i):endInd)));
            
            % if the max velo is maintained, and we didn't exceed the
            % midind
            if narrowFlag
                if endMaxVelo == currentMaxVelo && crossingStructIndexInBound(i) >= midInd ...
                        && jointAngleData.time(endInd) - jointAngleData.time(midInd) > 0.75
                    endIndToKeep = crossingStructIndexInBound(i);
                end
            else
                endIndToKeep = midInd;
            end
         
        end
        
        closeVal1 = jointAngleData.time(startIndToKeep);
        closeInd1 = startIndToKeep;
        closeVal2 = jointAngleData.time(endIndToKeep);
        closeInd2 = endIndToKeep;
        
%         [closeVal1, closeInd] = findClosestValue(startIndToKeep, jointAngleData.time);
%         [closeVal2, closeInd] = findClosestValue(endIndToKeep, jointAngleData.time);
    else
        % if no values are selected...
        [closeVal1, closeInd1] = findClosestValue(currTimeStart, jointAngleData.time); % insert a value so it's not zero. worst comes to worst, it's a value already used
        [closeVal2, closeInd2] = findClosestValue(currTimeEnd, jointAngleData.time); % insert a value so it's not zero. worst comes to worst, it's a value already used

    end
end

% function  [closeVal, closeInd] = ...
%             findZVCMid(currTimeStart, currTimeEnd, crossingStruct, jointAngleData, toleranceGap)
%     
%     midVal = 0.5*(currTimeStart + currTimeEnd);
%     lowerBound = currTimeStart < crossingStruct.Time;
%     upperBound = crossingStruct.Time < currTimeEnd;
%     inBound = and(lowerBound, upperBound);
% 
%     if sum(inBound)
%         % if there is a value
%         [midVal, midInd] = findClosestValue(midVal, crossingStruct.Time(inBound)); % pick the one closest to the middle
%         [closeVal, closeInd] = findClosestValue(midVal, jointAngleData.time);
%     else
%         % if no values are selected...
%         [closeVal, closeInd] = findClosestValue(currTimeStart, jointAngleData.time); % insert a value so it's not zero. worst comes to worst, it's a value already used
%     end
% end