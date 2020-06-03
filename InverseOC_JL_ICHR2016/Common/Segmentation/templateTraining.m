function [trainingData, dimReduct, classifier, classifierTrainingLabel, score, timeTaken] = ...
    templateTraining(dataSelect, dimReductSelect, classifierSelect, aggregatorSelect)
    % this function loads the appropriate data and trains the template.
    % it returns a dim reduction class and a classification class that is
    % used for classification purposes
    
    % -------------------- Data Selection --------------------
%     if 0
        trainingData = dataSelectProcessing(dataSelect);
%         save('trainingDataSave', 'trainingData');
%     else
%         load trainingDataSave
%     end
    
    % -------------------- Dim Reduction --------------------
    tic; % start timing after data has been loaded
    
    % train the dim reduct algorithm
    dimReduct = dimReductProcessing(dimReductSelect, trainingData);

    % run it on the training data for the classifier
%      trainingData.trainingDataDR = dimReduct.apply(trainingData.trainingData); % applying to the whole thing at once
       
    % apply PCA to each motionset 
    reducedData = [];
    for ind_motionSet = 1:length(trainingData.trainingIndMapping)
%         [closeVal, endingInd] = findClosestValue(trainingData.dataDivide(ind_motionSet), trainingData.trainingInd, 'below');

        if ind_motionSet == 1
            startingInd = 1;
            endingInd = trainingData.trainingIndMapping(ind_motionSet);
        else
            startingInd = trainingData.trainingIndMapping(ind_motionSet-1)+1;
            endingInd = trainingData.trainingIndMapping(ind_motionSet);
        end

        motionSetData = trainingData.trainingData(startingInd:endingInd, :);
        if strcmp(dimReductSelect.name, 'ksPCA')        
            motionSetDataDR = dimReduct.apply(motionSetData);
        else            
            motionSetDataDR = dimReduct.apply(motionSetData);
        end
%         motionSetDataDRZeroMean = motionSetDataDR;
        reducedData = [reducedData; motionSetDataDR];
    end
    
    trainingData.trainingDataDR = reducedData;
    
    % -------------------- Classifier --------------------
    [baseClassifier, baseClassifierLabel] = ...
        classifierProcessing(classifierSelect, aggregatorSelect, trainingData, dimReduct);
    
    groundTruthLabel = trainingData.trainingLabel;
    
    % pass the FSM from the training data to the classifier
    if isa(baseClassifier, 'SVMIncClassifier') 
        baseClassifier.adaptiveComponent(trainingData);
        
    elseif isa(baseClassifier, 'clusterClassifier') 
        baseClassifier.adaptiveComponent(trainingData);
        baseClassifierLabel = baseClassifier.classify(trainingData.trainingData, trainingData, 0);
        
    elseif isa(baseClassifier, 'PolyPCAClassifier')
        [baseClassifierLabel, groundTruthLabel] = baseClassifier.classify(trainingData.trainingData, trainingData);
        
    elseif isa(baseClassifier, 'SVMThreeClassClassifier')
        baseClassifier.adaptiveComponent(trainingData);
        baseClassifierLabel = baseClassifier.classify(trainingData.trainingDataDR);
        
        % for three_class, swap out training/testing_class3 with the
        % original label, since that's our actual assessment criteria
        trainingData.trainingLabel = trainingData.trainingLabel_class3;
    end
    
    % -------------------- Aggregator --------------------        
    [classifier, classifierTrainingLabel] = ...
        aggregatorProcessing(aggregatorSelect, baseClassifier, baseClassifierLabel, ...
            trainingData.trainingDataDR, trainingData.trainingLabel);
    
    timeTaken = toc;
        
    fprintf('Training Metrics (Tested only on points used for training) \n');
    score = accuracyAssess(classifierTrainingLabel, groundTruthLabel, trainingData.trainingName, ...
        trainingData.exerciseStruct.exerciseType, trainingData.exerciseStruct.exerciseType); % last one was 'dataSelect.patientSelection.fullExercistList' for a while
    printAccuracyToConsole(score);
    
    if 0
        figure;
        p1points = reducedData(trainingData.trainingLabel == 1, :);
        p1ds = randperm(size(p1points, 1), 2000);
        
        p0points = reducedData(trainingData.trainingLabel == 0, :);
        p0ds = randperm(size(p0points, 1), 2000);
        
        scatter3(p1points(p1ds, 1), p1points(p1ds, 2), p1points(p1ds, 3), 'bo')
        hold on; title('spca x1000')
        scatter3(p0points(p0ds, 1), p0points(p0ds, 2), p0points(p0ds, 3), 'rx')
    end
    
%     figure
%     clf
%     plot(classifierTrainingLabel,'bx');
%     hold on
%     plot(trainingData.trainingLabel-0.1,'rx');
%     plot(trainingData.trainingData(:, [1 4]) - 0.3);
%     title('red = testingLabel, blue = classifierLabel');
%     hold off
end