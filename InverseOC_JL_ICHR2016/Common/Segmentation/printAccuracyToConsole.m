function printAccuracyToConsole(score)
    % take in the results produced by accuracyAssess and print them to
    % screen
    fprintf('     TP: %u/%u (%f) \n', score.segmentTruePos, score.segmentTruePosTotal, score.segmentTruePosPercentage);
    fprintf('     TN: %u/%u (%f) \n', score.segmentTrueNeg, score.segmentTrueNegTotal, score.segmentTrueNegPercentage);
    fprintf('     FN: %u\n', score.segmentFalsePos);
    fprintf('     FP: %u\n', score.segmentFalseNeg);
    fprintf('     F1: %f (S) vs %f (C) \n', score.fScore_Seg, score.fScore_Class);
    fprintf('    Acc: %f (n) vs %f (b) \n', score.accuracy, score.bAcc);
    fprintf('    MCC: %f\n', score.MCC);
    fprintf('    AUC: %f\n', score.AUC);
end