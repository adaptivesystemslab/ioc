function [testingData, testingLabel, testingInd, testingName, testingTime, testingSegTime, testingIndMapping] = ...
    samplingForTestingData(obj, permissibleIndStruct, rawSegTime, rawPxLength, rawSegTestingInclude)

    sampleRules.px = [];
    sampleRules.tightSampling = 0;
    sampleRules.cropFirstLastPoints = 0;

    permissibleIndStruct.veloThreshold = 1:length(obj.label);
    permissibleIndStruct.downSample = 1:length(obj.label);
    permissibleIndStruct.SEP = 1:length(obj.label);
    permissibleIndStruct.segInclude = find(rawSegTestingInclude == 1); % exclude points that are on the don't include list
    
    if ~isfield(obj.settings, 'mode') || strcmpi(obj.settings.mode, 'Testing')
        % set the obs (though, if it's in 'both' mode,
        returnLimit = length(obj.label);

        [testingDataAll, testingLabelAll, testingIndAll, testingNameAll] = ...
            pullDataForTrainingPatients(obj, sampleRules, [], obj.settings.dimStrack, permissibleIndStruct, returnLimit);

        testingData = testingDataAll;
        testingLabel = testingLabelAll;
        testingInd = testingIndAll;
        testingName = testingNameAll;
        testingTime = obj.time(testingIndAll);
        testingSegTime = rawSegTime;
        testingIndMapping = testingIndAll;
    else
        testingData = [];
        testingLabel = [];
        testingInd = [];
        testingName = {};
        testingTime = [];
        testingSegTime = [];
        testingIndMapping = [];
    end
end