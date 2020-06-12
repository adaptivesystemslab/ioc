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