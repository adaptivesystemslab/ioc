function [fScore_Seg, fScore_Class] = fScoreCalculate(TP, TN, FP, FN, f1_beta)
    % calculate the F-score based on the inserted metrics
    % this function has been superceded by scoreCalculate, 2015-04-23
    
    correct = TP;
	fScore_Seg = (1+f1_beta^2)*correct ./ ((1+f1_beta^2)*correct + (f1_beta^2)*FN + FP);
    
    correct = TP + TN;
    fScore_Class = (1+f1_beta^2)*correct ./ ((1+f1_beta^2)*correct + (f1_beta^2)*FN + FP);
end