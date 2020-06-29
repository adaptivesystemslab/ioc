function [trainingData, trainingLabel, trainingInd, trainingName, trainingTime, trainingSegTime, trainingIndMapping] = ...
    samplingForTrainingData(obj, permissibleIndStruct, rawSegTime, rawPxLength, rawSegTrainingInclude)

    sampleRules.px = [];
    sampleRules.tightSampling = 0;
    sampleRules.cropFirstLastPoints = 0;
    
    if ~isfield(obj.settings, 'mode') || strcmpi(obj.settings.mode, 'Training')
        % pull out all the 'segment' and 'not-segment' entries
        returnLimit = floor(length(find(obj.label == obj.label_segment)) * obj.allowThresholdMultiplier);

        % make sure the threshold of total number of points are not
        % exceeded
        if returnLimit > obj.returnThreshold
            returnLimit = obj.returnThreshold;
        end

        sampleRules.tightSampling = 0;
        sampleRules.cropFirstLastPoints = 1;
        permissibleIndStruct.SEP = permissibleIndStruct.ind_allowSEPSegment;
        permissibleIndStruct.segInclude = find(rawSegTrainingInclude == obj.label_segment);
        [trainingDataSegment, trainingLabelSegment, trainingIndSegment, trainingNameSegment] = ...
            pullDataForTrainingPatients(obj, sampleRules, obj.label_segment, obj.settings.dimStrack, permissibleIndStruct, returnLimit);

        sampleRules.tightSampling = 0;
        sampleRules.cropFirstLastPoints = 0;
        permissibleIndStruct.SEP = 1:length(obj.label);
        permissibleIndStruct.segInclude = find(rawSegTrainingInclude ~= -1);
        [trainingDataUnknown, trainingLabelUnknown, trainingIndUnknown, trainingNameUnknown] = ...
            pullDataForTrainingPatients(obj, sampleRules, obj.label_unknown, obj.settings.dimStrack, permissibleIndStruct, returnLimit);

        
        sampleRules.tightSampling = obj.p0_guassianDownsampling;
        sampleRules.cropFirstLastPoints = 0;
        permissibleIndStruct.SEP = permissibleIndStruct.ind_allowSEPNotSegment;
        permissibleIndStruct.segInclude = find(rawSegTrainingInclude == obj.label_notSegment);
        [trainingDataNotSegment, trainingLabelNotSegment, trainingIndNotSegment, trainingNameNotSegment] = pullDataForTrainingPatients(obj, sampleRules, obj.label_notSegment, obj.settings.dimStrack, permissibleIndStruct, returnLimit);

        trainingDataTemp = [trainingDataSegment; trainingDataNotSegment; trainingDataUnknown];
        trainingLabelTemp = [trainingLabelSegment; trainingLabelNotSegment; trainingLabelUnknown];
        trainingIndTemp = [trainingIndSegment; trainingIndNotSegment; trainingIndUnknown];
        trainingNameTemp = [trainingNameSegment; trainingNameNotSegment; trainingNameUnknown];
        %                 trainingTime = obj.data(trainingIndTemp);

        [trainingIndSortedVal, trainingIndSortedInd] = sort(trainingIndTemp);

        for ind = 1:length(rawPxLength)
            [closeVal(ind), closeInd(ind)] = findClosestValue(rawPxLength(ind), trainingIndSortedVal, 'below');
        end

        trainingData = trainingDataTemp(trainingIndSortedInd, :);
        trainingLabel = trainingLabelTemp(trainingIndSortedInd);
        trainingInd = trainingIndSortedVal;
        trainingName = trainingNameTemp(trainingIndSortedInd);
        trainingTime = obj.time(trainingIndSortedVal);
        trainingSegTime = rawSegTime;
        trainingIndMapping = closeInd;
    else
        trainingData = [];
        trainingLabel = [];
        trainingName = {};
        trainingInd = [];
        trainingTime = [];
        trainingSegTime = [];
        trainingIndMapping = [];
    end
end