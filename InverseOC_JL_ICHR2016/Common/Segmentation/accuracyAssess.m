function score = accuracyAssess(classifierLabel, testingLabel, testingName, motionsToAssess, fullMotionList)
    % assess the accuracy of the labels reported
    
    % constants
    f1_beta = 1;
    longWay = 1; 
    
    % create a confusion matrix between the classes avaliable 
    possibleLabels = [0 1 2];
    settings.label_notSegment = 0;
    settings.label_segment = 1;
    settings.label_unknown = 2;
    
    notSegInd = find(possibleLabels == settings.label_notSegment);
    segInd = find(possibleLabels == settings.label_segment);
    unknownInd = find(possibleLabels == settings.label_unknown);
    
%     confMtxTemplate = zeros(length(possibleLabels)); % ground truth is y axis, classifier is x axis, motion is z axis
    
    % add in the 'non-motion' into all assessment
    if exist('motionsToAssess', 'var')
        motionsToAssess = ['-' motionsToAssess];
        
        if exist('fullMotionList', 'var')
            fullMotionList = ['-' fullMotionList];
        else
            fullMotionList = motionsToAssess;
        end
    end
    
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
            confMtxName(length(possibleLabels), length(possibleLabels), nameArrayInd) = 0; % create the array
        else
            nameArray = {};
            nameArrayInd = 0;
            tpArray = [];
            fnArray = [];
            fpArray = [];
            tnArray = [];
            confMtxName(length(possibleLabels), length(possibleLabels), 1) = 0; % create the array
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
                    confMtxName(2, 2, currNameInd) = 0; % create the array
                else
                    currNameInd = find(strcmpi(nameArray, currCurrName) == 1);
                end
                
                currClassifierLabel = classifierLabel(ind_name); % classifier label
                currTestingLabel = testingLabel(ind_name); % ground truth

                % ground truth is y axis, classifier is x axis, motion is z axis
                currClassLabelPos = find(currClassifierLabel == possibleLabels);
                currTestLabelPos = find(currTestingLabel == possibleLabels);
                confMtxName(currTestLabelPos, currClassLabelPos, currNameInd) = ...
                    confMtxName(currTestLabelPos, currClassLabelPos, currNameInd) + 1;
            end
        end
        
        [nameArraySort, nameArrayInd] = sort(nameArray);
        
        score.breakdownNameArray = nameArray(nameArrayInd)';
        score.breakdownTP = reshape(confMtxName(segInd, segInd, nameArrayInd), length(nameArray), 1, 1);
        score.breakdownTN = reshape(confMtxName(notSegInd, notSegInd, nameArrayInd), length(nameArray), 1, 1);
        score.breakdownFN = reshape(confMtxName(segInd, notSegInd, nameArrayInd), length(nameArray), 1, 1);
        score.breakdownFP = reshape(confMtxName(notSegInd, segInd, nameArrayInd), length(nameArray), 1, 1);
        
        % if the matrix is larger, we need a more comprehensive way to do
        % this
        
        
        breakdownTPSum = sum(score.breakdownTP);
        breakdownTNSum = sum(score.breakdownTN);
        breakdownFNSum = sum(score.breakdownFN);
        breakdownFPSum = sum(score.breakdownFP);

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
        confMtxAssessment = zeros(length(possibleLabels), length(possibleLabels), length(classifierLabel));
        
        for ind_assessment = 1:length(testingLabel)
            currClassifierLabel = classifierLabel(ind_assessment); % classifier label
            currTestingLabel = testingLabel(ind_assessment); % ground truth
            
            currClassLabelPos = find(currClassifierLabel == possibleLabels);
            currTestLabelPos = find(currTestingLabel == possibleLabels);
            confMtxAssessment(currTestLabelPos, currClassLabelPos, ind_assessment) = 1;
        end
        
        truePos_array = reshape(confMtxAssessment(segInd, segInd, :), length(testingLabel), 1, 1);
        trueNeg_array = reshape(confMtxAssessment(notSegInd, notSegInd, :), length(testingLabel), 1, 1);
        falseNeg_array = reshape(confMtxAssessment(segInd, notSegInd, :), length(testingLabel), 1, 1);
        falsePos_array = reshape(confMtxAssessment(notSegInd, segInd, :), length(testingLabel), 1, 1);
        
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
    
% % %     % indicies for checking
% % %     incorrect_segmentInd = find(falseNeg_array == 1); % 'segment', incorrect
% % %     incorrect_nonSegmentInd = find(falsePos_array == 1); % 'non-segment', incorrect
% % %     
% % %     numeric_array = 1:length(testingLabel);
% % %     ack1 = numeric_array(truePos_array == 1);
% % %     ack2 = find(classifierLabel == 1);
% % %     correct_segmentInd = intersect(ack1, ack2); % 'segment', correct
% % %     
% % %     ack1 = numeric_array(truePos_array == 1);
% % %     ack2 = find(classifierLabel == 0);
% % %     correct_nonSegmentInd = intersect(ack1, ack2); % 'non-segment', correct
    
    score.segmentTruePos = truePos;
    score.segmentTruePosTotal = truePosTotal;
    score.segmentTruePosPercentage = truePos/truePosTotal;
    score.segmentTrueNeg = trueNeg;
    score.segmentTrueNegTotal = trueNegTotal;
    score.segmentTrueNegPercentage = trueNeg/trueNegTotal;
    score.segmentFalsePos = falsePos;
    score.segmentFalseNeg = falseNeg;
    
    % F1 score
    score = scoreCalculate(score, f1_beta);
%     [fScore_Seg, fScore_Class, accuracy, MCC, bAcc] = fScoreCalculate(truePos, trueNeg, falsePos, falseNeg, f1_beta);
    
    % calculate AUC
    if ~isempty(testingLabel) && length(unique(testingLabel)) > 1
        [X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(testingLabel,classifierLabel,1);
  %     plot(X, Y); hold on; plot([0 1], [0 1], 'r--');
    else
        AUC = 0;
    end

    score.AUC = AUC;
    score.confMtx = sum(confMtxAssessment, 3);
    score.indAssessed = indAssessed;
    
    if exist('testingName', 'var')
        % need to wait for the main function to finish it's calculations
        % before we can calculate the rest, for the F1Total value
        header = {'Motion', 'TP', 'TN', 'FN', 'FP', 'F1_Seg', 'F1_Class'};
        F1Sum = {'Total (Double counts datapoints that has two labels)', ...
            breakdownTPSum, breakdownTNSum, breakdownFNSum, breakdownFPSum, fScoreBreakdown_Seg, fScoreBreakdown_Class};
        F1Total = {'Total (Does not double count)', truePos, trueNeg, falseNeg, falsePos, score.fScore_Seg, score.fScore_Class};
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