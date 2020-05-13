function [classifierTestingLabel, classifierTestingProb, testingData, groundTruthLabel, groundTruthNames] ...
    = adaptiveProcessing(dataSelect, dimReductSelect, classifierSelect, aggregatorSelect, ...
    trainingData, motionSetData, dimReduct, classifier, testingData)

    % processing the adaptive components of the dimreduct/classifier
    % updates, if applicable
    
    downsampleUpperBound = 5000; % so n p1 points, and n p0 points
    
    % either 'apca' or 'akspca'
    if strcmp(dimReductSelect.settingName(1), 'a')
        % keep the original training data on top of the new data
        fprintf('Retraining aPCA/Classifier...\n');
        
        % apply downsample and save the training data
        if size(trainingData.trainingData, 1)/2 > downsampleUpperBound
            % randomly select some points for downsampling
            p1points = find(trainingData.trainingLabel == 1); % don't want to assume the array is sorted
            p0points = find(trainingData.trainingLabel == 0);
            
            selectInd = randperm(length(p1points), downsampleUpperBound); % first half is p1 points
            selectInd2 = randperm(length(p0points), downsampleUpperBound);
            
            X_tr = trainingData.trainingData([p1points(selectInd) p0points(selectInd2)], :);
            Y_tr = trainingData.trainingLabel([p1points(selectInd) p0points(selectInd2)], :);
        else
            X_tr = trainingData.trainingData;
            Y_tr = trainingData.trainingLabel;
        end
        
        newData.trainingLabel = [Y_tr; testingData.testingLabel];
        newData.trainingData = [X_tr; testingData.testingData];
        dimReduct = dimReductProcessing(dimReductSelect, newData);
        
        switch classifierSelect.name
            case 'SVMInc'
                % don't forcibly update the SVMInc since it might be set to
                % retain data from previous runs
                
                % instead, update its rotation matrix
                classifier.updateDimReductClass(dimReduct);
                
                % pass the FSM from the training data to the classifier
                classifier.adaptiveComponent(trainingData);
            
            otherwise
                % apply the new PCA rot to the full original training
                trainingData.trainingDataDR = dimReduct.apply(trainingData.trainingData);
                        
                % now also need to retrain the classifier with new PCA 
                % (but not the obs data)
                [classifier, ~] = ...
                    classifierProcessing(classifierSelect, aggregatorSelect, trainingData);
        end
    else
        % not adapting the dim reduction, so don't need to retrain the
        % classifier or update any dimreduct class in inc svm
    end
    
    % weither or not the dim reduct is updated, we want to apply the
    % transformation
    switch dimReductSelect.name
        case 'ksPCA'
            motionSetDataDR = dimReduct.apply(motionSetData);
            
        otherwise
            motionSetDataDR = dimReduct.apply(motionSetData);
    end
    testingData.testingDataDR = motionSetDataDR;
    
    charString = [testingData.subjectData{1}.datasetName ...
        '_Subj' num2str(testingData.subjectData{1}.subjectNumber) ...
        '_Sess' num2str(testingData.subjectData{1}.sessionNumber) ...
        '_' testingData.subjectData{1}.exerciseName];
    
    testingData.identString = charString;
    
    groundTruthLabel = testingData.testingLabel;
    groundTruthNames = testingData.testingName;
    
    % run the classifier
    if isa(classifier, 'SVMIncClassifier') || isa(classifier, 'clusterClassifier') 
        classifier.initializeAdaptiveElements;
        [classifierTestingLabel, classifierTestingProb] = classifier.classify(testingData.testingData, testingData, 1);
        
    elseif isa(classifier, 'PolyPCAClassifier')
        [classifierTestingLabel, groundTruthLabel] = classifier.classify(testingData.testingData, testingData, 1);
        groundTruthNames = testingData.testingName(1:length(groundTruthLabel));
        
    elseif isa(classifier, 'SVMThreeClassClassifier')
        classifierTestingLabel = classifier.classify(testingData.testingDataDR);
        
        % for three_class, swap out training/testing_class3 with the
        % original label, since that's our actual assessment criteria
        testingData.testingLabel = testingData.testingLabel_class3;
    else
        classifierTestingLabel = classifier.classify(testingData.testingDataDR);
        classifierTestingProb = 0.5 * ones(size(classifierTestingLabel));
    end
end