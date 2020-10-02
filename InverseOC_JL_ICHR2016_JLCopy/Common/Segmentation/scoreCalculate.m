function overallScore = scoreCalculate(instanceScore, f1_beta)
    % calculate various scoring metrics for assessment purposes

    instanceScore = formatInstanceScore(instanceScore);
    
    overallScore = ...
        fScoreCalculate_instance(instanceScore.segmentTruePos, instanceScore.segmentTrueNeg, ...
        instanceScore.segmentFalsePos, instanceScore.segmentFalseNeg, f1_beta);

    % merge the original struct with the new one generated from _instance
    names = [fieldnames(instanceScore); fieldnames(overallScore)];
    overallScore = cell2struct([struct2cell(instanceScore); struct2cell(overallScore)], names, 1);
    
    % sum all the arrays together
    overallScore.segmentTruePosTotal = sum(instanceScore.segmentTruePosTotal);
    overallScore.segmentTruePos = sum(instanceScore.segmentTruePos);
    overallScore.segmentTruePosPercentage = overallScore.segmentTruePos/overallScore.segmentTruePosTotal;
    overallScore.segmentTrueNegTotal = sum(instanceScore.segmentTrueNegTotal);
    overallScore.segmentTrueNeg = sum(instanceScore.segmentTrueNeg);
    overallScore.segmentTrueNegPercentage = overallScore.segmentTrueNeg/overallScore.segmentTrueNegTotal;
    overallScore.segmentFalsePos = sum(instanceScore.segmentFalsePos);
    overallScore.segmentFalseNeg = sum(instanceScore.segmentFalseNeg);
end

function instanceScore = formatInstanceScore(instanceScore)
    if ~isfield(instanceScore, 'segmentTruePos') && isfield(instanceScore, 'TP')
        instanceScore.segmentTruePos = instanceScore.TP;
        instanceScore.segmentTrueNeg = instanceScore.TN;
        instanceScore.segmentFalsePos = instanceScore.FP;
        instanceScore.segmentFalseNeg = instanceScore.FN;
        instanceScore.segmentTruePosTotal = instanceScore.TPTotal;
        instanceScore.segmentTrueNegTotal = instanceScore.TNTotal;
    end
end

function [fScore_Seg, fScore_Class, accuracy, MCC, bAcc] = fScoreCalculate_overall(TP, TN, FP, FN, f1_beta)
    % calculate the F-score based on the inserted metrics. this is the
    % original metric calculation method 
    correct = TP;
	fScore_Seg = (1+f1_beta^2)*correct ./ ((1+f1_beta^2)*correct + (f1_beta^2)*FN + FP);
    
    correct = TP + TN;
    fScore_Class = (1+f1_beta^2)*correct ./ ((1+f1_beta^2)*correct + (f1_beta^2)*FN + FP);
    
    % calculate the accuracy
    accuracy = (TP + TN) ./ (TP + TN + FP + FN);
    
    % calculate the MCC
    num = (TP .* TN) - (FP .* FN);
    den = (TP + FP) .* (TP + FN) .* (TN + FP) .* (TN + FN);
    
    den(den == 0) = 1;

    MCC = num ./ sqrt(den);
    
    % balanced accuracy
    if TP == 0 && FN == 0
        bAcc_1 = 0;
    else
        bAcc_1 = TP/(TP+FN);
    end
    
    if TN == 0 && FP == 0
        bAcc_2 = 0;
    else
        bAcc_2 = TN/(TN+FP);
    end
    
    bAcc = (1/2)*(bAcc_1 + bAcc_2);
end

function assessment = fScoreCalculate_instance(TP, TN, FP, FN, f1_beta)
    % calculate the F-score based on the inserted metrics, by weighing the
    % individual values
    
    TPTotal = TP+FN;
    TNTotal = TN+FP;
    TotalPoints = TPTotal+TNTotal;
    OVerallTotal = sum(TotalPoints);
    weights = TotalPoints/OVerallTotal; % multiply everything by this
    
    correct = TP;
	fScore_Seg = (1+f1_beta^2)*correct ./ ((1+f1_beta^2)*correct + (f1_beta^2)*FN + FP);
    
    correct = TP + TN;
    fScore_Class = (1+f1_beta^2)*correct ./ ((1+f1_beta^2)*correct + (f1_beta^2)*FN + FP);
    
    % calculate the accuracy
    accuracy = (TP + TN) ./ (TP + TN + FP + FN);
    
    % calculate the MCC
    num = (TP .* TN) - (FP .* FN);
    den = (TP + FP) .* (TP + FN) .* (TN + FP) .* (TN + FN);
    
    den(den == 0) = 1;

    MCC = num ./ sqrt(den);
    
    % balanced accuracy
    bAcc_1 = TP./(TP+FN);
    bAcc_2 = TN./(TN+FP);
    
    bAcc_1(isnan(bAcc_1)) = 0;
    bAcc_2(isnan(bAcc_2)) = 0;
    
    bAcc = (1/2)*(bAcc_1 + bAcc_2);
    
    if isnan(fScore_Seg)
        assessment.fScore_Seg = 0;
        assessment.fScore_Seg_std = 0;
        assessment.fScore_Class = 0;
        assessment.fScore_Class_std = 0;
        assessment.accuracy = 0;
        assessment.accuracy_std = 0;
        assessment.MCC = 0;
        assessment.MCC_std = 0;
        assessment.bAcc = 0;
        assessment.bAcc_std = 0;
    else
        [assessment.fScore_Seg, assessment.fScore_Seg_std]     = sumStd(fScore_Seg, weights);
        [assessment.fScore_Class, assessment.fScore_Class_std] = sumStd(fScore_Class, weights);
        [assessment.accuracy, assessment.accuracy_std]         = sumStd(accuracy, weights);
        [assessment.MCC, assessment.MCC_std]                   = sumStd(MCC, weights);
        [assessment.bAcc, assessment.bAcc_std]                 = sumStd(bAcc, weights);
    end

end

function [sumVal, stddevVal] = sumStd(baseVal, weights)
    weightedValues = baseVal .* weights;
    
    sumVal = sum(weightedValues);    
    stddevVal = sqrt(var(baseVal, weights));
end