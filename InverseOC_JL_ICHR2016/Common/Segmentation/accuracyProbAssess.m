function score = accuracyProbAssess(classifierLabel, classifierProb, testingLabel, testingName, motionsToAssess, fullMotionList)
    % assess the accuracy of the labels reported
    
    % constants
    f1_beta = 1;
    longWay = 1; 
    
    % add in the 'non-motion' into all assessment
    if exist('motionsToAssess', 'var')
        motionsToAssess = ['-' motionsToAssess];
        
        if exist('fullMotionList', 'var')
            fullMotionList = ['-' fullMotionList];
        else
            fullMotionList = motionsToAssess;
        end
    end
    
    classifierProb = classifierProb - 0.5; % anything above zero = class 1
    
    % calculate the breakdown if the proper data is passed in. This section
    % double counts datapoints that have two labels. the subsequent section
    % single counts datapoints regardless of the label count
    if exist('testingName', 'var')
        indAssessed = zeros(size(classifierLabel)); % track which data points are being including in the scoring
        
        if exist('fullMotionList', 'var')
            % if there is a predefined motion list we want
            nameArray = fullMotionList;
            nameArrayInd = length(nameArray);
            tpArray(nameArrayInd) = 0;
            fnArray(nameArrayInd) = 0;
            fpArray(nameArrayInd) = 0;
            tnArray(nameArrayInd) = 0;
            
            tpArray_prob(nameArrayInd, 2) = 0; % high, low prob
            fnArray_prob(nameArrayInd, 2) = 0;
            fpArray_prob(nameArrayInd, 2) = 0;
            tnArray_prob(nameArrayInd, 2) = 0;
        else
            nameArray = {};
            nameArrayInd = 0;
            tpArray = [];
            fnArray = [];
            fpArray = [];
            tnArray = [];
            
            tpArray_prob = []; % high, low prob
            fnArray_prob = []; % high, low prob
            fpArray_prob = []; % high, low prob
            tnArray_prob = []; % high, low prob
        end
        
        for ind_name = 1:length(testingName)
            % insert the name into the array
            currName = testingName{ind_name};
            if isempty(currName)
                currName = {'-'}; % label the no-segments as something
            elseif ischar(currName)
                temp = currName;
                currName = {}; % define currName as a cell first, if it's a single entry
                currName{1} = temp;
            end
            
            % remove repeats, though am counting all the
%             other items twice, so perhaps we won't remove
%             currName = unique(currName); 
            
            for ind_currName = 1:length(currName)
                currCurrName = currName{ind_currName};
                
                if exist('motionsToAssess', 'var') && ...
                        sum(strcmpi(currCurrName, motionsToAssess)) == 0
                    % motion is not in our 'approved' list of motions to
                    % check for
                    continue
                end
                
                % passed all the tests, start tallying
                indAssessed(ind_name) = 1;

                if sum(strcmpi(nameArray, currCurrName)) == 0
                    nameArrayInd = nameArrayInd + 1;
                    nameArray{nameArrayInd} = currCurrName;
                    currNameInd = nameArrayInd;
                    tpArray(currNameInd) = 0;
                    fnArray(currNameInd) = 0;
                    fpArray(currNameInd) = 0;
                    tnArray(currNameInd) = 0;
                    
                    tpArray_prob(currNameInd, 2) = 0; % high, low prob
                    fnArray_prob(currNameInd, 2) = 0;
                    fpArray_prob(currNameInd, 2) = 0;
                    tnArray_prob(currNameInd, 2) = 0;
                else
                    currNameInd = find(strcmpi(nameArray, currCurrName) == 1);
                end
                
                currClassifierLabel = classifierLabel(ind_name);
                currClassifierProb = abs(classifierProb(ind_name));
                currTestingLabel = testingLabel(ind_name);
                
                classifierThreshold = 0.3;

                if currClassifierLabel == 1 && currTestingLabel == 1
                    % true positive
                    tpArray(currNameInd) = tpArray(currNameInd) + 1;
                    
                    if currClassifierProb > classifierThreshold
                        tpArray_prob(currNameInd, 1) = tpArray_prob(currNameInd, 1) + 1;
                    else
                        tpArray_prob(currNameInd, 2) = tpArray_prob(currNameInd, 2) + 1;
                    end
                    
                elseif currClassifierLabel == 1 && currTestingLabel == 0
                    % false positive
                    fpArray(currNameInd) = fpArray(currNameInd) + 1;
                    
                    if currClassifierProb > classifierThreshold
                        fpArray_prob(currNameInd, 1) = fpArray_prob(currNameInd, 1) + 1;
                    else
                        fpArray_prob(currNameInd, 2) = fpArray_prob(currNameInd, 2) + 1;
                    end
                    
                elseif currClassifierLabel == 0 && currTestingLabel == 1
                    % false negative
                    fnArray(currNameInd) = fnArray(currNameInd) + 1;
                    
                    if currClassifierProb > classifierThreshold
                        fnArray_prob(currNameInd, 1) = fnArray_prob(currNameInd, 1) + 1;
                    else
                        fnArray_prob(currNameInd, 2) = fnArray_prob(currNameInd, 2) + 1;
                    end
                else
                    tnArray(currNameInd) = tnArray(currNameInd) + 1;
                    
                    if currClassifierProb > classifierThreshold
                        tnArray_prob(currNameInd, 1) = tnArray_prob(currNameInd, 1) + 1;
                    else
                        tnArray_prob(currNameInd, 2) = tnArray_prob(currNameInd, 2) + 1;
                    end
                end
            end
        end
        
        [nameArraySort, nameArrayInd] = sort(nameArray);
        
        score.breakdownNameArray = nameArray(nameArrayInd)';
        score.breakdownTP = tpArray(nameArrayInd)';
        score.breakdownTN = tnArray(nameArrayInd)';
        score.breakdownFN = fnArray(nameArrayInd)';
        score.breakdownFP = fpArray(nameArrayInd)';
        
        score.breakdownTP_higProb = tpArray_prob(nameArrayInd, 1)';
        score.breakdownTP_lowProb = tpArray_prob(nameArrayInd, 2)';
        score.breakdownTN_higProb = tnArray_prob(nameArrayInd, 1)';
        score.breakdownTN_lowProb = tnArray_prob(nameArrayInd, 2)';
        score.breakdownFN_higProb = fnArray_prob(nameArrayInd, 1)';
        score.breakdownFN_lowProb = fnArray_prob(nameArrayInd, 2)';
        score.breakdownFP_higProb = fpArray_prob(nameArrayInd, 1)';
        score.breakdownFP_lowProb = fpArray_prob(nameArrayInd, 2)';
        
        breakdownTPSum = sum(score.breakdownTP);
        breakdownTNSum = sum(score.breakdownTN);
        breakdownFNSum = sum(score.breakdownFN);
        breakdownFPSum = sum(score.breakdownFP);
        
        breakdownTPSum_higProb = sum(score.breakdownTP_higProb);
        breakdownTPSum_lowProb = sum(score.breakdownTP_lowProb);
        breakdownTNSum_higProb = sum(score.breakdownTN_higProb);
        breakdownTNSum_lowProb = sum(score.breakdownTN_lowProb);
        breakdownFNSum_higProb = sum(score.breakdownFN_higProb);
        breakdownFNSum_lowProb = sum(score.breakdownFN_lowProb);
        breakdownFPSum_higProb = sum(score.breakdownFP_higProb);
        breakdownFPSum_lowProb = sum(score.breakdownFP_lowProb);

        [fScoreIndividual_Seg, fScoreIndividual_Class] = fScoreCalculate(score.breakdownTP, score.breakdownTN, score.breakdownFP, score.breakdownFN, f1_beta); %((1+f1_beta^2)*breakdownT./((1+f1_beta^2)*breakdownT+(f1_beta^2)*score.breakdownFN+score.breakdownFP));
        [fScoreBreakdown_Seg, fScoreBreakdown_Class] = fScoreCalculate(breakdownTPSum, breakdownTNSum, breakdownFPSum, breakdownFNSum, f1_beta);
    else
        % breakdown data is not passed in. perform assessment based on just
        % the classifier and ground truth arrays
        indAssessed = ones(size(classifierLabel));
    end
    
    % crop down the classifier and testing label accordingly. we are either
    % assessing the accuracy of the whole classifierLabel array, or just a
    % subset, if we've only been calculating a subset from above. Either
    % case, the above 'testingName' section double counts all datapoints
    % that has more than one label (weighing on the individual labels),
    % while this section only counts each data point as a single point,
    % regardless of how many labels there are (ie similar to how a
    % classifier should be assessed)
    classifierLabel = classifierLabel(indAssessed == 1);
    testingLabel = testingLabel(indAssessed == 1);
    
    if ~longWay
        % pre-setup matrices        
        truePos_array = and(classifierLabel, testingLabel);
        falseNeg_array = and(xor(classifierLabel, testingLabel), classifierLabel);
        falsePos_array = and(xor(classifierLabel, testingLabel), testingLabel);
    else
        scoreArray = zeros(length(classifierLabel), 4);
        
        for ind_assessment = 1:length(testingLabel)
            if classifierLabel(ind_assessment) == 1 && testingLabel(ind_assessment) == 1
                % true positive
                scoreArray(ind_assessment, 1) = 1;
            elseif classifierLabel(ind_assessment) == 1 && testingLabel(ind_assessment) == 0
                % false positive
                scoreArray(ind_assessment, 2) = 1;
            elseif classifierLabel(ind_assessment) == 0 && testingLabel(ind_assessment) == 1
                % false negative
                scoreArray(ind_assessment, 3) = 1;
            elseif classifierLabel(ind_assessment) == 0 && testingLabel(ind_assessment) == 0
                % true negative
                scoreArray(ind_assessment, 4) = 1;
            else
                fprintf('Unknown label, not categorized. \n');
            end
        end
        
        truePos_array = scoreArray(:, 1); %
        trueNeg_array = scoreArray(:, 4);
        falseNeg_array = scoreArray(:, 3);
        falsePos_array = scoreArray(:, 2);
        
        correct_array = truePos_array + trueNeg_array;
        incorrect_array = falseNeg_array + falsePos_array;
    end
    
    truePosTotal = sum(testingLabel == 1);
    trueNegTotal = sum(testingLabel == 0);
    
    % tally
    truePos = sum(truePos_array);
    trueNeg = sum(trueNeg_array);
    falsePos = sum(falsePos_array);
    falseNeg = sum(falseNeg_array);
    
    % indicies for checking
    incorrect_segmentInd = find(falseNeg_array == 1); % 'segment', incorrect
    incorrect_nonSegmentInd = find(falsePos_array == 1); % 'non-segment', incorrect
    
    numeric_array = 1:length(testingLabel);
    ack1 = numeric_array(truePos_array == 1);
    ack2 = find(classifierLabel == 1);
    correct_segmentInd = intersect(ack1, ack2); % 'segment', correct
    
    ack1 = numeric_array(truePos_array == 1);
    ack2 = find(classifierLabel == 0);
    correct_nonSegmentInd = intersect(ack1, ack2); % 'non-segment', correct
    
    % F1 score
    [fScore_Seg, fScore_Class, accuracy, MCC, bAcc] = fScoreCalculate(truePos, trueNeg, falsePos, falseNeg, f1_beta);
    
    % calculate AUC
    if ~isempty(testingLabel) && length(unique(testingLabel)) > 1
        [X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(testingLabel,classifierLabel,1);
  %     plot(X, Y); hold on; plot([0 1], [0 1], 'r--');
    else
        AUC = 0;
    end
    
    score.indAssessed = indAssessed;
    score.segmentTruePos = truePos;
    score.segmentTruePosTotal = truePosTotal;
    score.segmentTruePosPercentage = truePos/truePosTotal;
    score.segmentTrueNeg = trueNeg;
    score.segmentTrueNegTotal = trueNegTotal;
    score.segmentTrueNegPercentage = trueNeg/trueNegTotal;
    score.segmentFalsePos = falsePos;
    score.segmentFalseNeg = falseNeg;
    score.fScore_Seg = fScore_Seg;
    score.fScore_Class = fScore_Class;
    score.accuracy = accuracy; 
    score.MCC = MCC;
    score.AUC = AUC; 
    score.bAcc = bAcc;
    score.breakdownTPSum_higProb = breakdownTPSum_higProb;
    score.breakdownTPSum_lowProb = breakdownTPSum_lowProb;
    score.breakdownTNSum_higProb = breakdownTNSum_higProb;
    score.breakdownTNSum_lowProb = breakdownTNSum_lowProb;
    score.breakdownFNSum_higProb = breakdownFNSum_higProb;
    score.breakdownFNSum_lowProb = breakdownFNSum_lowProb;
    score.breakdownFPSum_higProb = breakdownFPSum_higProb;
    score.breakdownFPSum_lowProb = breakdownFPSum_lowProb;
    
    if exist('testingName', 'var')
        % need to wait for the main function to finish it's calculations
        % before we can calculate the rest, for the F1Total value
        header = {'Motion', 'TP', 'TN', 'FN', 'FP', 'F1_Seg', 'F1_Class'};
        F1Sum = {'Total (Double counts datapoints that has two labels)', ...
            breakdownTPSum, breakdownTNSum, breakdownFNSum, breakdownFPSum, fScoreBreakdown_Seg, fScoreBreakdown_Class};
        F1Total = {'Total (Does not double count)', truePos, trueNeg, falseNeg, falsePos, fScore_Seg, fScore_Class};
        combinedArray = [score.breakdownTP score.breakdownTN ...
            score.breakdownFN score.breakdownFP ...
            fScoreIndividual_Seg fScoreIndividual_Class];
        totalCombined = [header; score.breakdownNameArray num2cell(combinedArray); ...
            F1Sum; F1Total];
        score.breakdownTotalCombined = totalCombined;
    end
    
    score.correct_array = correct_array;
    score.incorrect_array = incorrect_array;
%     score.correct_segmentInd = correct_segmentInd;
%     score.correct_nonSegmentInd = correct_nonSegmentInd;
%     score.incorrect_segmentInd = incorrect_segmentInd;
%     score.incorrect_nonSegmentInd = incorrect_nonSegmentInd;
end