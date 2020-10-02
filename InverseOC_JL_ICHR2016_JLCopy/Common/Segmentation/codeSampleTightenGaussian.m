function trainingInd = codeSampleTightenGaussian(trainingInd, labelTemp)
    % from intersectInd, we can find all the points that has a gap in it.
    gapToleranceThreshold = 0.1; % pull this much points from around the gap
    gapDownsampleThreshold = 0.3; % downsample the rest of the data

    diffTrainingArray = [1; diff(trainingInd)];
    gapTrainingArrayMid = find(diffTrainingArray > 1);

    gapTrainingArray = [1; gapTrainingArrayMid; length(trainingInd)];
    
%     if trainingInd(1) == 1
%         gapTrainingArray = [1; gapTrainingArray];
%     end

%     if trainingInd(end) == length(labelTemp)
%         gapTrainingArray(end+1) = length(labelTemp);
%     end

    trainingPrepArray = [];

    for ind_gap = 1:length(gapTrainingArray)-1
        % pull out the values close to the gaps
        startGap = gapTrainingArray(ind_gap);
        endGap = gapTrainingArray(ind_gap+1)-1;

        if endGap == startGap  || endGap+1 == startGap
            % no change distance between the gaps...how did this happen?
            continue
        end
        
        % pull out these cases now
        gapLength = endGap - startGap;
        sampleArray = 1:gapLength;
        gaussArray = gauss(0, 10*gapLength, sampleArray'); % runtime
%       gaussArray = gauss(0, 3*gapLength, sampleArray'); % plotting
        normalizedGaussArray = gaussArray/gaussArray(1);
        
        % generate the threshold pattern for the left and right gap
        thresholdLeft = rand(size(normalizedGaussArray));
        passVarLeft = normalizedGaussArray > thresholdLeft;
        
        thresholdRight = rand(size(normalizedGaussArray));
        passVarRight = flipud(normalizedGaussArray > thresholdRight);
        
        % combine the two
        passVar = or(passVarLeft, passVarRight);
        gapInd = startGap:endGap-1;
        gapKeep = gapInd(passVar)';
%         gapKeep = trainingInd(gapInd(passVar));

        trainingPrepArray = [trainingPrepArray; gapKeep;];
        
%         plot(passVarLeft, 'x')
%         hold on
%         plot(passVarRight, 'rx')
%         plot(passVar, 'go')
    end

    trainingInd = sort(trainingInd(trainingPrepArray));
end