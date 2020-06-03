function [scoreTemporal, h] = temporalClusterLocalization(saveStruct, dataSelect, exerciseSubsetList)

    if ~exist('exerciseSubsetList', 'var') || isempty(exerciseSubsetList)
        exerciseSubsetInd = 1:length(saveStruct.segTimeStartGroundTruth);
    else
        exerciseSubsetInd = [];
        for ind_segInclude = 1:length(saveStruct.segName)
            if sum(strcmpi(exerciseSubsetList, saveStruct.segName(ind_segInclude)))
                exerciseSubsetInd(end+1) = ind_segInclude;
            end
        end
    end

    saveStruct.segTimeStartGroundTruth = saveStruct.segTimeStartGroundTruth(exerciseSubsetInd);
    saveStruct.segTimeEndGroundTruth = saveStruct.segTimeEndGroundTruth(exerciseSubsetInd);
    saveStruct.segName = saveStruct.segName(exerciseSubsetInd);
    
    [manualSegStruct, algSegStruct, h] = convertClustersToSegmentPoints(saveStruct, ...
        dataSelect.temporalTolerance.tolerance, dataSelect.settings.segmentPointWindow, dataSelect.plotFig);

    % remove black list items

    % separate into classify and generalize for the comparison

    % select only
    metricTemporal = manualSegmentationComparison(manualSegStruct, algSegStruct, dataSelect.temporalTolerance);
    scoreTemporal.segmentTruePosTotal = metricTemporal.totalManualCount;
    scoreTemporal.segmentTruePos = metricTemporal.correct;
    scoreTemporal.segmentTrueNegTotal = 0;
    scoreTemporal.segmentTrueNeg = 0;
    scoreTemporal.segmentTrueNegPercentage = 0;
    scoreTemporal.segmentFalsePos = metricTemporal.falsePos;
    scoreTemporal.segmentFalseNeg = metricTemporal.falseNeg;
    scoreTemporal.segmentTruePosPercentage = scoreTemporal.segmentTruePos / scoreTemporal.segmentTruePosTotal;
    scoreTemporal.correctArray = metricTemporal.correctArray;
    scoreTemporal.falseNegArray = metricTemporal.falseNegArray;
    [scoreTemporal.fScore_Seg, scoreTemporal.fScore_Class] = fScoreCalculate(scoreTemporal.segmentTruePos, ...
        0, scoreTemporal.segmentFalsePos, scoreTemporal.segmentFalseNeg, 1);
end