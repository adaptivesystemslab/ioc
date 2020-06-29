function metric = manualSegmentationComparison(manualSegStruct, algSegStruct, sysparam, timeTolArray)
    % performs the algorithmic segmentation accurac check by iterating
    % through the manual segments and see if algorithmic segments fall
    % within the error bounds
    
    % current potiential limitations
    % - if small segment window is declared in the error tolerance of the
    % first segment point, it would consume that particular manual segment.
    % then, if a correct segment was also found within the that, it would
    % be declared incorrect. 
    
    if ~exist('sysparam', 'var')
        timeTolerance = (1/33.3333) * 8;
        timeOffset = -0;
        
        % so add the offset directly to the bounds, ie currSegTime = currSegTime + timeOffset;
        offsetDir = 'shift';
        
        %or 'expand', which apples it unevenly to the start and end
        % ie currSegTime(:, 1) - timeOffset; currSegTime(:, 2) + timeOffset;
%         offsetDir = 'expand';
    else
        timeTolerance = sysparam.tolerance;
        timeOffset = sysparam.offset;
        offsetDir = sysparam.offsetDir; 
    end
    
    constToIgnore = -1; % once a value is spent, it should be ignored
    smallPosFloat = -0.01; % if currSegTime has a value that's less than 0, it is set to this
    
    if exist('timeTolArray', 'var')
        timeTolInner = timeTolArray(1);
        timeTolOuter = timeTolArray(2);
    else

    % 1.95N on one, 0.05N on other
%     timeTolInner = 0.05;
%     timeTolOuter = 1.95;
    
    % 1N on one, 1N on other
    timeTolInner = 1;
    timeTolOuter = 1;
    
    % 1.50N on one, 0.5N on other
%     timeTolInner = 0.5;
%     timeTolOuter = 1.5;
    
    % 1.95N on one, 1N on other
%     timeTolInner = 1;
%     timeTolOuter = 1.95;
    
    % 1.95N on one, 1N on other
%     timeTolInner = 1.95;
%     timeTolOuter = 1.95;
    end

    [manualSegTime, mapping] = sort(manualSegStruct.manualSegTime);
    manualSegLabel = manualSegStruct.label(mapping);
    manualSegLabelMap = sort(unique(manualSegStruct.label));
    
    algSegTime = algSegStruct.segmentTime;
    algSegLabel = algSegStruct.segmentLabel;
    
    algSegTimeOrig = algSegTime;
    
    counterManualCountTotal = zeros(1, length(manualSegLabelMap));
    counterAlgCountTotal = zeros(1, length(manualSegLabelMap));
    
    counterCorrect = zeros(1, length(manualSegLabelMap));
    counterFalseNeg = zeros(1, length(manualSegLabelMap));
    counterFalseNeg_tooFar = zeros(1, length(manualSegLabelMap));
    counterFalseNeg_missingSegment = zeros(1, length(manualSegLabelMap));
    counterFalsePos = zeros(1, length(manualSegLabelMap));
    counterFalsePos_extraSegment = zeros(1, length(manualSegLabelMap));
    
    arrayCorrect = cell(1, length(manualSegLabelMap));
    arrayFalseNeg = cell(1, length(manualSegLabelMap));
    arrayFalseNegBase = cell(1, length(manualSegLabelMap));
    arrayFalsePos = cell(1, length(manualSegLabelMap));

    for i = 1:length(manualSegLabelMap)
        activeManualLabel = manualSegLabelMap{i};
        
        activeCounterCorrect = 0;
        activeCounterFalseNeg = 0;
        activeCounterFalseNeg_tooFar = 0;
        activeCounterFalseNeg_missingSegment = 0;        
        activeCounterFalsePos = 0;
        activeCounterFalsePos_extraSegment = 0;
        
        activeArrayCorrect = []; % correct entries
        activeArrayFalseNeg = []; % entries that should be there but off
        activeArrayFalseNegBase = []; % where it should've been
        activeArrayFalsePos = []; % extra segments that should not exist
        activeArrayManCorrect = []; % correct entries, but from man side
        
        % pull all the segments with this label (manual)
        activeManualLabelArrayInd = find(strcmpi(manualSegLabel, activeManualLabel));
        activeManualSegTime = manualSegTime(activeManualLabelArrayInd);
        activeManualCountTotal = length(activeManualLabelArrayInd);
        
        % now pull all the segments reported from the algorithm
        activeAlgLabelArrayInd = find(strcmpi(algSegLabel, activeManualLabel));
        activeAlgSegTime = algSegTime(activeAlgLabelArrayInd);
        activeAlgCountTotal = length(activeAlgLabelArrayInd);
        
        % iterating through the manual seg points and see if there are
        % corresponding alg seg points 
        for j = 1:2:length(activeManualSegTime)             
            if j == 1
                % first run, so no previous segment is declared
                prevSegTime = [smallPosFloat; smallPosFloat; smallPosFloat];
            else
                prevSegTime = currSegTime(:, 2);
            end
            
            currSegTime = [activeManualSegTime(j)-timeTolerance*timeTolOuter activeManualSegTime(j+1)-timeTolerance*timeTolInner; ...
                           activeManualSegTime(j)                            activeManualSegTime(j+1); ...
                           activeManualSegTime(j)+timeTolerance*timeTolInner activeManualSegTime(j+1)+timeTolerance*timeTolOuter];

           if strcmp(offsetDir, 'expand')
               currSegTime(:, 1) = currSegTime(:, 1) + timeOffset(1);
               currSegTime(:, 2) = currSegTime(:, 2) + timeOffset(2);
           else
               currSegTime = currSegTime + timeOffset;
           end
           
            zeroChecker = find(currSegTime <= 0);
            if ~isempty(zeroChecker)
                currSegTime(zeroChecker) = smallPosFloat;
            end
                       
            % a flag that a pending FN (missing entry) is coming up. if 
            % nothing clears it by the end of the cycle, declare an FN
            declareFalseNegative = 0;
            declareFalseNegativeArray = []; % the ind1 and ind2 for it
                       
            % the gap before the segment - if both points are before the
            % segment, then we have a false positive, an extra segment that
            % is declared that should not have been declared
            % CHECK - upper bound of prevSegTime to lower bound of currSegTime(j)
            preSegInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'preSegInd');
            
            while ~isempty(preSegInd)
                ind1 = preSegInd(1);
                ind2 = preSegInd(1)+1;
                segPt1 = activeAlgSegTime(ind1); % this one is less than lowerbound of currSegTime(j)
                segPt2 = activeAlgSegTime(ind2); % not sure which this one is yet
                
                if segPt2 < currSegTime(3, 1) 
                    % also less than lowerbound of currSegTime(j+1)
                    % 2x FALSE POSITIVE (extra segment)
                    activeArrayFalsePos = [activeArrayFalsePos segPt1 segPt2];
                    
                    activeCounterFalsePos = activeCounterFalsePos + 2;
                    activeCounterFalsePos_extraSegment = activeCounterFalsePos_extraSegment + 2;

                elseif segPt2 < currSegTime(1, 2)
                    % in the gap between the segment points...
                    % 2x FALSE NEGATIVE (both segment point wrong)
                    activeArrayFalseNeg = [activeArrayFalseNeg segPt1 segPt2];
                    activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, :)];
                    
                    activeCounterFalseNeg = activeCounterFalseNeg + 2;
                    activeCounterFalseNeg_tooFar = activeCounterFalseNeg_tooFar + 2;
                    
                elseif segPt2 < currSegTime(3, 2)
                    % if there is a single unconsumed segment point, it 
                    % means that it's overlapping in the preseg region and 
                    % in the seg if it is less than the lower bound of the 
                    % end segment point, it is totally in the middle, and 
                    % should be marked as a false negative
                    
                    % so between lower and upper bound of (j+1)
                    % 1x FALSE NEGATIVE (first segment point wrong)
                    % 1X CORRECT(second segment point correct)
                    activeArrayFalseNeg = [activeArrayFalseNeg segPt1];
                    activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, 1)];
                    activeArrayCorrect = [activeArrayCorrect segPt2];
                    activeArrayManCorrect = [activeArrayManCorrect currSegTime(2, 2)];
                    
                    activeCounterCorrect = activeCounterCorrect + 1;
                    activeCounterFalseNeg = activeCounterFalseNeg + 1;
                    activeCounterFalseNeg_tooFar = activeCounterFalseNeg_tooFar + 1;
                    
                else
                    % beyond upperbound of (j+1)
                    % 2x FALSE NEGATIVE (both segment point wrong)
                    activeArrayFalseNeg = [activeArrayFalseNeg segPt1 segPt2];
                    activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, :)];
                    
                    activeCounterFalseNeg = activeCounterFalseNeg + 2;
                    activeCounterFalseNeg_tooFar = activeCounterFalseNeg_tooFar + 2;
                end
                
                activeAlgSegTime(ind1) = constToIgnore;
                activeAlgSegTime(ind2) = constToIgnore;
                
                preSegInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'preSegInd');
            end
            
            % by now, all the presegment points should be consumed, with
            % exception to a currSegTime(i) being too early but has its
            % currSegTime(i+1) within the second error bound
            
            % check for first segment point error region 
            segInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'segInd');
            
            % prevent double counting of a false neg from preSeg section
            intersectArray = intersect(activeArrayFalseNegBase, currSegTime(2, :));
            
            if length(segInd) == 0 && isempty(intersectArray)
                % nothing detected in the region
                % 2x FALSE NEGATIVE (missing segment)
%                 activeArrayFalseNeg = [activeArrayFalseNeg segPt1 segPt2];

                % don't want to declare a FN just yet, in case there's some
                % delayed points later
                declareFalseNegative = 1;
                declareFalseNegativeArray = currSegTime(2, :); % the ind1 and ind2 for it
            
%                 activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, :)];
%                 
%                 activeCounterFalseNeg = activeCounterFalseNeg + 2;
%                 activeCounterFalseNeg_missingSegment = activeCounterFalseNeg_missingSegment + 2;
                
            elseif length(segInd) > 0
                % correctly has at least one segment point
                ind1 = segInd(1); % good point
                ind2 = segInd(1)+1; % not sure yet
                segPt1 = activeAlgSegTime(ind1); % this one is less than lowerbound of currSegTime(j)
                segPt2 = activeAlgSegTime(ind2); % not sure which this one is yet
               
                % 1x CORRECT(first segment point correct)
                activeArrayCorrect = [activeArrayCorrect segPt1];
                activeArrayManCorrect = [activeArrayManCorrect currSegTime(2, 1)];
                activeCounterCorrect = activeCounterCorrect + 1;
                
                if segPt2 < currSegTime(1, 2) 
                    % second point is too close to first                   
                    % 1x FALSE NEGATIVE (second segment point wrong)                    
                    activeArrayFalseNeg = [activeArrayFalseNeg segPt2];
                    activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, 2)];
                    
                    activeCounterFalseNeg = activeCounterFalseNeg + 1;
                    activeCounterFalseNeg_tooFar = activeCounterFalseNeg_tooFar + 1;
                    
                elseif segPt2 < currSegTime(3, 2)
                    % second point is within second bound
                    % 1x CORRECT (second point correct)
                    
                    activeArrayCorrect = [activeArrayCorrect segPt2];
                    activeArrayManCorrect = [activeArrayManCorrect currSegTime(2, 2)];
                    activeCounterCorrect = activeCounterCorrect + 1;
                else
                    % and too far away
                    % 1x FALSE NEGATIVE (second segment point wrong)
                    activeArrayFalseNeg = [activeArrayFalseNeg segPt2];
                    activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, 2)];
                    
                    activeCounterFalseNeg = activeCounterFalseNeg + 1;
                    activeCounterFalseNeg_tooFar = activeCounterFalseNeg_tooFar + 1;
                end
                
                activeAlgSegTime(ind1) = constToIgnore;
                activeAlgSegTime(ind2) = constToIgnore;
                
                % now check the entire segment bracket. there should not be
                % anything else here,, since we already picked up a segment
                % from the first section...but if there is, then they're bad
                segBadInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'segBadInd');
                while ~isempty(segBadInd)
                    ind1 = segBadInd(1); % bad point, excess segment
                    ind2 = segBadInd(1)+1; % also bad. bad bad
                    segPt1 = activeAlgSegTime(ind1);
                    segPt2 = activeAlgSegTime(ind2);
                    
                    % 2x FALSE POSITIVE (extra segments)
                    activeArrayFalsePos = [activeArrayFalsePos segPt1 segPt2];
                    
                    activeCounterFalsePos = activeCounterFalsePos + 2;
                    activeCounterFalsePos_extraSegment = activeCounterFalsePos_extraSegment + 2;
                    
                    activeAlgSegTime(ind1) = constToIgnore;
                    activeAlgSegTime(ind2) = constToIgnore;
                    
                    segBadInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'segBadInd');
                end
            end
            
            % checking the space between the first and second segment, to
            % the end of the second segment. there shouldn't be anything
            % here. there could be a case of too large first segment but
            % correct second segment
            % TODO1: need to decide if should be currSegTime(1, 2) or currSegTime(3, 2)...
            segBadInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'segBadIndBetween');
            
%             activeArrayManCorrect = [activeArrayManCorrect currSegTime(2, 2)];
            while ~isempty(segBadInd)
                ind1 = segBadInd(1); % bad point, excess segment
                ind2 = segBadInd(1)+1; % also bad. bad bad
                segPt1 = activeAlgSegTime(ind1);
                segPt2 = activeAlgSegTime(ind2);
                
                % once in here, need to make sure that the (first) point
                % we're looking at isn't actually part of the next segment
                if j+2 < length(activeManualSegTime) % so not at the end yet
                    if activeManualSegTime(j+2) - timeTolerance < segPt1
                        % it is within the next bracket...so save it for
                        % the next bracket
                        break
                    end
                end
                
                % prevent double counting of a correct neg from seg section
                intersectArray = intersect(activeArrayManCorrect, currSegTime(2, :));
                
                if segPt2 < currSegTime(1, 2)
                    % too close in the middle
                    % 2x FALSE POSITIVE (extra segments)
                    activeArrayFalsePos = [activeArrayFalsePos segPt1 segPt2];
                    
                    activeCounterFalsePos = activeCounterFalsePos + 2;
                    activeCounterFalsePos_extraSegment = activeCounterFalsePos_extraSegment + 2;
                    
                elseif segPt2 < currSegTime(3, 2) && isempty(intersectArray)
                    % 1x FALSE NEGATIVE (first point too far)
                    % 1x CORRECT (second segment correct. current segment not included in correct yet)
                    activeArrayCorrect = [activeArrayCorrect segPt2];
                    activeArrayManCorrect = [activeArrayManCorrect currSegTime(2, 2)];
                    activeArrayFalseNeg = [activeArrayFalseNeg segPt1];
                    activeArrayFalseNegBase = [activeArrayFalseNegBase currSegTime(2, 1)];
                    
                    activeCounterCorrect = activeCounterCorrect + 1;
                    activeCounterFalseNeg = activeCounterFalseNeg + 1;
                    activeCounterFalseNeg_tooFar = activeCounterFalseNeg_tooFar + 1;
                    
                    % earlier, since the first point isn't in segInd
                    % region, a FN flag would've been raised. clear the
                    % flag now, since we've gotten a match
                    declareFalseNegative = 0;
                    declareFalseNegativeArray = [];
                    
                else
                    % 2x FALSE POSITIVE (extra segments)
                    activeArrayFalsePos = [activeArrayFalsePos segPt1 segPt2];
                    
                    activeCounterFalsePos = activeCounterFalsePos + 2;
                    activeCounterFalsePos_extraSegment = activeCounterFalsePos_extraSegment + 2;
                end
                
                activeAlgSegTime(ind1) = constToIgnore;
                activeAlgSegTime(ind2) = constToIgnore;
                
                segBadInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'segBadIndBetween');
            end
            
            if declareFalseNegative
                % an FN flag was raised. consume the data
                activeArrayFalseNegBase = [activeArrayFalseNegBase declareFalseNegativeArray];
                
                activeCounterFalseNeg = activeCounterFalseNeg + 2;
                activeCounterFalseNeg_missingSegment = activeCounterFalseNeg_missingSegment + 2;
                
                % just in case something bad happens, make sure it's cleared                
                declareFalseNegative = 0;
                declareFalseNegativeArray = [];
            end
        end
        
        % checking after the last declared segment point 
        postSegInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'postSegInd');
        
        while ~isempty(postSegInd)
            ind1 = postSegInd(1); % bad point, excess segment
            ind2 = postSegInd(1)+1; % also bad. bad bad
            segPt1 = activeAlgSegTime(ind1);
            segPt2 = activeAlgSegTime(ind2);
            
            % 2x FALSE POSITIVE (extra segments)
            activeArrayFalsePos = [activeArrayFalsePos segPt1 segPt2];
            
            activeCounterFalsePos = activeCounterFalsePos + 2;
            activeCounterFalsePos_extraSegment = activeCounterFalsePos_extraSegment + 2;
            
            activeAlgSegTime(ind1) = constToIgnore;
            activeAlgSegTime(ind2) = constToIgnore;
            
            postSegInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, 'postSegInd');
        end
        
        counterManualCountTotal(i) = activeManualCountTotal;
        counterAlgCountTotal(i) = activeAlgCountTotal;
        
        counterCorrect(i) = activeCounterCorrect;
        
        counterFalseNeg(i) = activeCounterFalseNeg;
        counterFalseNeg_tooFar(i) = activeCounterFalseNeg_tooFar;
        counterFalseNeg_missingSegment(i) = activeCounterFalseNeg_missingSegment;
        
        counterFalsePos(i) = activeCounterFalsePos;        
        counterFalsePos_extraSegment(i) = activeCounterFalsePos_extraSegment;

        arrayCorrect{i} = activeArrayCorrect;
        arrayFalseNeg{i} = activeArrayFalseNeg;
        arrayFalseNegBase{i} = activeArrayFalseNegBase;
        arrayFalsePos{i} = activeArrayFalsePos;
    end
    
    % calculate the number of templates that appear in alg, but is not part
    % of the manual segment at all (ie identified incorrect templates)
    oddOneOutAlgIds = setxor(manualSegStruct.id, algSegStruct.segmentId);
    falsePos_oddOneOutCounter = 0;
    for i = 1:length(oddOneOutAlgIds)
        falsePos_oddOneOutCounter = falsePos_oddOneOutCounter + length(find(oddOneOutAlgIds(i) == algSegStruct.segmentId));
    end
    
    metric.order = manualSegLabelMap;
    metric.timeTolerance = timeTolerance;
    metric.timeTolOffset = [timeTolInner timeTolOuter];
    
    metric.totalAlgCount = counterAlgCountTotal;
    metric.totalManualCount = counterManualCountTotal;
    
    metric.correct = counterCorrect;
    metric.falsePos = counterFalsePos;
    metric.falsePos_extra = counterFalsePos_extraSegment;
    metric.falseNeg = counterFalseNeg;
    metric.falseNeg_tooFar = counterFalseNeg_tooFar;
    metric.falseNeg_missing = counterFalseNeg_missingSegment;
    
    metric.totalManualCountSum = sum(counterManualCountTotal);
    metric.totalAlgCountSum = sum(counterAlgCountTotal);
    
    metric.correctSum = sum(counterCorrect);
    metric.falsePosSum = sum(counterFalsePos) + falsePos_oddOneOutCounter;
    metric.falsePos_extraSum = sum(counterFalsePos_extraSegment);
    metric.falsePos_nonTemplateSet = falsePos_oddOneOutCounter;
    metric.falseNegSum = sum(counterFalseNeg);
    metric.falseNeg_tooFarSum = sum(counterFalseNeg_tooFar);
    metric.falseNeg_missingSum = sum(counterFalseNeg_missingSegment);
    
    num = metric.correctSum + metric.falseNeg_tooFarSum;
    denom = metric.correctSum + metric.falsePosSum + metric.falseNegSum;
    metric.accuracy_correctonly = metric.correctSum / denom;
    metric.accuracy_FNincluded = num / denom;

    metric.correctArray = arrayCorrect;
    metric.falsePosArray = arrayFalsePos;
    metric.falseNegArray = arrayFalseNeg;
    metric.falseNegBaseArray = arrayFalseNegBase;
    
    breakdownMtx = [metric.totalManualCount;
        metric.correct;
        metric.falsePos;
        metric.falseNeg;
        metric.totalAlgCount];
    
    summaryMtx = [metric.totalManualCountSum; 
        metric.correctSum;
        metric.falsePosSum;
        metric.falseNegSum;
        metric.falsePos_nonTemplateSet;
        metric.falseNeg_tooFarSum;
        metric.falseNeg_missingSum];
    
    metric.summaryMtx = sum(summaryMtx, 2);
    metric.breakdownMtx = breakdownMtx;
end

function [counterCorrectRevised, counterFalsePosRevised, uniqueArrayCorrect, uniqueArrayFalsePosXOR] = overlapCorrection(arrayCorrect, arrayFalsePos, counterCorrect, counterFalsePos)

    uniqueArrayCorrectCombine = [];
    for i = 1:length(arrayCorrect)
        uniqueArrayCorrect{i} = sort(unique(arrayCorrect{i})); % unique entries + zeros are removed
        
        if ~isempty(uniqueArrayCorrect{i}) && uniqueArrayCorrect{i}(1) == 0
            % if there are zeros, it is removed here
            uniqueArrayCorrect{i} = uniqueArrayCorrect{i}(2:end);
        end
        
        % updating the new counter
        counterCorrectRevised(i) = length(uniqueArrayCorrect{i});
        
        % combining it into a large array for later comparisons
        uniqueArrayCorrectCombine = [uniqueArrayCorrectCombine; uniqueArrayCorrect{i}];
    end
    uniqueArrayCorrectCombine = sort(uniqueArrayCorrectCombine);
    
    % sorting out the false positive array as well, and removing zeros    
    uniqueArrayFalsePos = sort(unique(arrayFalsePos));

    if ~isempty(uniqueArrayFalsePos) && uniqueArrayFalsePos(1) == 0
        uniqueArrayFalsePos = uniqueArrayFalsePos(2:end);
    end
    
    % any entries that is in the false pos and also in correct may be
    % double counted, so this will remove such double counted entries
    arrayIntersect = intersect(uniqueArrayCorrectCombine, uniqueArrayFalsePos);
    uniqueArrayFalsePosXOR = setxor(uniqueArrayFalsePos, arrayIntersect);
    counterFalsePosRevised = length(uniqueArrayFalsePosXOR);
end

function segInd = segIndFinder(activeAlgSegTime, prevSegTime, currSegTime, indType)
    % making a function for this so each of the segment finding are done 
    % uniform, since there are a lot of different checks being done, a 
    % modification in one area may not be carried out properly
    
    switch indType
        case 'preSegInd'
            % the gap before the segment - if both points are before the
            % segment, then we have a false positive, an extra segment that
            % is declared that should not have been declared
            % CHECK - upper bound of prevSegTime to lower bound of currSegTime(j)
            % TODO1: prevSegTime(1) vs prevSegTime(3)? (check other TODO1)
%             segInd = find(prevSegTime(1) < activeAlgSegTime & activeAlgSegTime < currSegTime(1, 1));
            segInd = find(0 <= activeAlgSegTime & activeAlgSegTime < currSegTime(1, 1));
            
        case 'segInd'
            % check for first segment point error region
            segInd = find(currSegTime(1, 1) <= activeAlgSegTime & activeAlgSegTime <= currSegTime(3, 1));
            
        case 'segBadInd'
            % now check the entire segment bracket. there should not be
            % anything else here, since we already picked up a segment
            % from the first section...but if there is, then they're bad
            segInd = find(currSegTime(1, 1) <= activeAlgSegTime & activeAlgSegTime <= currSegTime(3, 1));
            
        case 'segBadIndBetween'
            % checking the space between the first and second segment pt, to
            % the end of the second segment. there shouldn't be anything
            % here. there could be a case of too large first segment but
            % correct second segment
            % TODO1: need to decide if should be currSegTime(1, 2) or currSegTime(3, 2)...
            segInd = find(currSegTime(3, 1) <= activeAlgSegTime & activeAlgSegTime <= currSegTime(1, 2)); 
            
        case 'postSegInd'
            % checking after the last declared segment point 
            segInd = find(activeAlgSegTime > currSegTime(1, 2));
    end
end	

function metric = manualSegmentationComparison_old2(manualSegStruct, algSegTime, algSegId, timeTolerance)
    % manualSegmentationComparison(algSegTime, 1, 'R', manualSegTime)
    % settings
%     mode = 1; % 1 = single (up/down together)
%     mode = 2; % 2 = double (up/down separate)
    
%     dataType = '3';
%     dataType = 'R';

%     timeOffsetThres = 0.3; % allow for this much error to pass metric, for start/end time
%     timeLengthThres = 0.3; % allow for this much error to pass metric, for segment length
    if ~exist('timeTolerance', 'var')
        timeTolerance = (1/60) * 12;
    end
    
    [manualSegTime, mapping] = sort(manualSegStruct.Time);
    manualSegLabel = manualSegStruct.Label(mapping);
    manualSegLabelMap = sort(unique(manualSegStruct.Label));
    
    algSegTimeOrig = algSegTime;
    
    counterCorrect = zeros(1, length(manualSegLabelMap));
    counterFalseNeg = zeros(1, length(manualSegLabelMap));
    counterFalsePos = 0;
    counterFalsePosBin = zeros(1, length(manualSegLabelMap));
    
    arrayCorrect = cell(1, length(manualSegLabelMap));
    arrayFalseNeg = cell(1, length(manualSegLabelMap));
    arrayFalsePos = [];
    arrayFalsePosBin = cell(1, length(manualSegLabelMap));

    prevSegTime = -timeTolerance;
    for i = 1:length(manualSegTime)        
        currSegTime = manualSegTime(i);
        currSegTimeLower = currSegTime - timeTolerance;
        currSegTimeHigher = currSegTime + timeTolerance;
        prevSegTimeHigher = prevSegTime + timeTolerance;
        
        % the gap before the segment point
        % - false positive (if there are seg points)
        % check between upper range of prev seg point and lower range of
        % current seg point
        preSegInd = find(prevSegTimeHigher < algSegTime & algSegTime < currSegTimeLower);
        
        if length(preSegInd) > 0
            arrayFalsePos = [arrayFalsePos; preSegInd'];
            counterFalsePos = counterFalsePos + length(preSegInd);
        end
        
        % segment point
        % - false negative (no seg points)
        % - correct (exactly one point)
        % - false positive (too many seg points)
        % check between the lower range of current seg point and upper
        % range of current seg point
        segInd = find(currSegTimeLower < algSegTime & algSegTime < currSegTimeHigher);
        segPointsInRange = length(segInd);
        currentManSegLabel = manualSegLabel{i};
        currentManSegInd = find(strcmpi(manualSegLabelMap, currentManSegLabel));
        
        if segPointsInRange == 0
            % no alg seg detected. false negative
            counterFalseNeg(currentManSegInd) = counterFalseNeg(currentManSegInd) + 1;
            arrayFalseNeg{currentManSegInd} = [arrayFalseNeg{currentManSegInd}; currSegTime];
        else
            % there are alg seg detected
            if strcmpi(currentManSegLabel, algSegId{segInd(1)})
                counterCorrect(currentManSegInd) = counterCorrect(currentManSegInd) + 1;
                
                
%                 if isempty(arrayCorrect)
%                     arrayCorrect = [arrayCorrect; segInd(1)];
%                 elseif arrayCorrect(end) == segInd(1)
%                     arrayCorrect = [arrayCorrect; segInd(2)];
%                     segInd(2) = 0;
%                 else
%                     arrayCorrect = [arrayCorrect; segInd(1)];
%                     segInd(1) = 0;
%                 end
                
                arrayCorrect{currentManSegInd} = [arrayCorrect{currentManSegInd}; segInd(1)];
                algSegTime(segInd(1)) = 0;
                

%                 [closeVal, closeInd] = findClosestValue(currSegTime, algSegTime(segInd));
%                 
%                 arrayCorrect = [arrayCorrect; segInd(closeInd)];
%                 segInd(closeInd) = 0;
            else
                counterFalsePos = counterFalsePos + 1;
                arrayFalsePos = [arrayFalsePos; segInd(1)];
            end
            
%             if segPointsInRange > 0 && segInd(1) ~= 0
            if segPointsInRange > 1
                counterFalsePos = counterFalsePos + (segPointsInRange - 1);
                arrayFalsePos = [arrayFalsePos; segInd(2:end)'];
            end
        end
        
        prevSegTime = manualSegTime(i);
    end
    
    % space after last segment point (false positives)
    LastSegTimeHigher = manualSegTime(end) + timeTolerance;
    postSegInd = find(algSegTime > LastSegTimeHigher);
    counterFalsePos = counterFalsePos + length(postSegInd);
    arrayFalsePos = [arrayFalsePos; postSegInd'];
    
    [counterCorrectRevised, counterFalsePosRevised, arrayCorrectRevised, arrayFalsePosRevised] = overlapCorrection(arrayCorrect, arrayFalsePos, counterCorrect, counterFalsePos);
     
    % sort out the remaining false positives
    for i = 1:length(arrayFalsePosRevised)
        ind = arrayFalsePosRevised(i);
        currSegTime = algSegTime(ind);
        currentManSegInd = find(strcmpi(manualSegLabelMap, algSegId{ind}));
        arrayFalsePosBin{currentManSegInd} = [arrayFalsePosBin{currentManSegInd}; currSegTime];
        counterFalsePosBin(currentManSegInd) = counterFalsePosBin(currentManSegInd) + 1;
    end
    
%     metric.segcountMan = manualSegCount;
%     metric.segcountAlg = algSegCount;
    metric.correct = counterCorrectRevised;
    metric.correctSum = sum(counterCorrectRevised);
    
    metric.falseNeg = counterFalseNeg;
    metric.falseNegSum = sum(counterFalseNeg);
    
    metric.falsePos = counterFalsePosBin;
    metric.falsePosSum = sum(counterFalsePosBin);
    
    metric.correctArray = arrayCorrectRevised;
    metric.falsePosArray = arrayFalsePosRevised;
    metric.falseNegArray = arrayFalsePosBin;
end

function metric = manualSegmentationComparison_old1(algSegTime, mode, dataType, manualSegTime)
    % manualSegmentationComparison(algSegTime, 1, 'R', manualSegTime)
    % settings
%     mode = 1; % 1 = single (up/down together)
%     mode = 2; % 2 = double (up/down separate)
    
%     dataType = '3';
%     dataType = 'R';

    timeOffsetThres = 0.3; % allow for this much error to pass metric, for start/end time
    timeLengthThres = 0.3; % allow for this much error to pass metric, for segment length

    if ~exist('manualSegTime', 'var')
        % load manual segmentation data
        basePath = 'C:\Users\jf2lin\Documents\MATLAB\APARS\Segmentation\ManualSegmentation\';
        baseFile = 'jlin_legExt_singleSegment.mat';
        jlinTemplate = load([basePath baseFile]);
        
        if dataType == '3'
            manualSegTime = jlinTemplate.segInfo3.t;
            manualSegCount = jlinTemplate.segInfo3.segCount;
        elseif dataType == 'R'
            manualSegTime = jlinTemplate.segInfoR.t;
            manualSegCount = jlinTemplate.segInfoR.segCount;
        end
        
        manualSegTimeOdd = manualSegTime(1:2:end);
        manualSegTimeEven = manualSegTime(2:2:end);
        manualTimeLength = manualSegTimeEven - manualSegTimeOdd;
    end
    
    switch mode
        case 1
            % up/down together
            segModifier = 2;
            
        case 2
            % up/down separate
            segModifier = 4;
    end
    
    % check for number of segmentation
    algSegCount = floor(length(algSegTime)/segModifier);
    countDelta = manualSegCount - algSegCount;

    % check for starting and ending positioning for segments
    algTimeMatch = length(algSegTime);
    algTimeLength = length(algSegTime/2);
    for i = 1:algSegCount
        oddInd = i*2 - 1;
        evenInd = i*2;
        
        [algTimeMatchStart(i), algTimeMatchStartInd(i)] = findClosestValue(algSegTime(oddInd), manualSegTimeOdd);
        [algTimeMatchEnd(i), algTimeMatchEndInd(i)] = findClosestValue(algSegTime(evenInd), manualSegTimeEven);
        algTimeLength(i) = algSegTime(evenInd) - algSegTime(oddInd);
        algTimeMatch(oddInd) = algTimeMatchStart(i);
        algTimeMatch(evenInd) = algTimeMatchEnd(i);
    end
    algTimeMatch = algTimeMatch';
    algTimeLength = algTimeLength';
    timeDelta = abs(algTimeMatch - algSegTime); 
    timeDeltaFind = find(timeDelta > timeOffsetThres);
    
    % check for length of segments
    for i = 1:algSegCount
        lengthDelta(i) = abs(algTimeLength(i) - manualTimeLength(algTimeMatchStartInd(i)));
    end
%     lengthDelta = abs(manualTimeLength - algTimeLength);
    lengthDeltaFind = find(lengthDelta > timeLengthThres);
    
    metric.segcountMan = manualSegCount;
    metric.segcountAlg = algSegCount;
    metric.timeDeltaFind = length(timeDeltaFind);
    metric.lengthDeltaFind = length(lengthDeltaFind);
end