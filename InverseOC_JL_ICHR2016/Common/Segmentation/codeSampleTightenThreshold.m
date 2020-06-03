function trainingInd = codeSampleTightenThreshold(trainingInd, labelTemp)
    % from intersectInd, we can find all the points that has a gap in it.
    gapToleranceThreshold = 0.1; % pull this much points from around the gap
    gapDownsampleThreshold = 0.3; % downsample the rest of the data

    diffTrainingArray = [1; diff(trainingInd)];
    gapTrainingArray = find(diffTrainingArray > 1);

    if trainingInd(1) == 1
        gapTrainingArray = [1; gapTrainingArray];
    end

    if trainingInd(end) == length(labelTemp)
        gapTrainingArray(end+1) = length(labelTemp);
    end

    trainingPrepArray = [];

    for ind_gap = 1:length(gapTrainingArray)-1
        % pull out the values close to the gaps
        startGap = gapTrainingArray(ind_gap);
        endGap = gapTrainingArray(ind_gap+1)-1;

        % pull out these cases now
        gapLength = endGap - startGap;
        gapSkim = floor(gapLength*gapToleranceThreshold);

        if startGap+gapSkim > length(trainingInd)
            continue
        end

        gapStartKeep = trainingInd(startGap:startGap+gapSkim);
        gapEndKeep = trainingInd(endGap-gapSkim:endGap);

        % downsample the remaining values
        gapDownsample = trainingInd(startGap+gapSkim+1:endGap-gapSkim-1);
        gapDownsampleAllInd = randperm(length(gapDownsample));
        gapDownsampleSelectInd = gapDownsampleAllInd(1:floor(length(gapDownsample)*gapDownsampleThreshold));
        gapDownsampleKeep = gapDownsample(gapDownsampleSelectInd);

        trainingPrepArray = [trainingPrepArray; gapStartKeep; gapDownsampleKeep; gapEndKeep];
    end

    trainingInd = sort(trainingPrepArray);
end