function runEval_single(dataSelect, dimReductSelect, classifierSelect, aggregatorSelect) 
    % take in a series of structures designed to denote the specific
    % combination of data to load
    
    if ~exist('dataSelect', 'var')
        % default settings
        dataSelect.name = 'Healthy1'; % Cosine, Healthy1
        
%         dimReductSelect.name = 'FDA'; % None, PCA, FDA
%         classifierSelect.name = 'RBF'; % kNN, QDA, RBF, SVM, NN
%         aggregatorSelect.name = 'None'; % None, Bagging, Boosting
        
        dimReductSelect.name = 'PCA'; % None, PCA, FDA
        classifierSelect.name = 'SVM'; % kNN, QDA, RBF, SVM, NN
        aggregatorSelect.name = 'Boosting'; % None, Bagging, Boosting
    end
    
    % -------------------- Data Selection --------------------
    dataSelect.settings.label_notSegment = 0;
    dataSelect.settings.label_segment = 1;
    dataSelect.settings.segmentPointWindow = 4; % how much non-segment points to convert to segment point
    dataSelect.settings.dimStrack = 9; % set to 0 if only want a single point per row
    
    data = dataSelectProcessing(dataSelect);
    
    % -------------------- Dim Reduction --------------------
    fprintf('Training Dim Reduction...\n');
    [dimReduct] = dimReductProcessing(dimReductSelect, data);
    
    data.trainingDataDR = dimReduct.apply(data.trainingData);
    data.testingDataDR = dimReduct.apply(data.testingData);
    
    % -------------------- Classifier and Aggregation--------------------
    fprintf('Training Classifier...\n');
    [classifier, classifierTrainingLabel] = ...
        classifierTraining(aggregatorSelect, classifierSelect, data.trainingDataDR, data.trainingLabel);
    
    classifierTestingLabel = classifier.classify(data.testingDataDR);
    
    % -------------------- Post-processing --------------------    
    
    figure
    hold on
    plot(data.testingLabel, 'rx');
    plot(classifierTestingLabel - 0.1, 'bo');    
    ylim([-0.5 1.5]);
    title('red = testing label, blue = classifier label');
    
    % perform accuracy analysis
    score = accuracyAssess(classifierTrainingLabel, data.trainingLabel);
    fprintf('Training Metrics\n');
    fprintf('Correct: %u/%u (%f) \n', score.segmentCorrect, score.segmentCorrectTotal, score.segmentCorrectPercentage);
    fprintf('     FN: %u\n', score.segmentFalsePos);
    fprintf('     FP: %u\n', score.segmentFalseNeg);
    
    score = accuracyAssess(classifierTestingLabel, data.testingLabel);
    fprintf('Testing Metrics\n');
    fprintf('Correct: %u/%u (%f)\n', score.segmentCorrect, score.segmentCorrectTotal, score.segmentCorrectPercentage);
    fprintf('     FN: %u\n', score.segmentFalseNeg);
    fprintf('     FP: %u\n', score.segmentFalsePos);
end

% % % function data = dataSelectProcessing(dataSelect)
% % %     switch dataSelect.name
% % %         case 'Cosine'
% % %             if ~isfield(dataSelect, 't')
% % %                 % generate default settings
% % %                 dataSelect.t = 0:pi/32:8*pi;
% % %                 dataSelect.amp = 1;
% % %                 dataSelect.phase = 0;
% % %             end              
% % %             
% % %             data = CosineData(dataSelect.settings, dataSelect.t, dataSelect.amp, dataSelect.phase);
% % %             
% % %             
% % %         case 'Healthy1'
% % %             if ~isfield(dataSelect, 'basePath')
% % %                 % generate default settings
% % %                 dataSelect.basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy1_2011-11\';
% % %                 dataSelect.patientSelection.pxCount = [3];
% % %                 dataSelect.phase = 0;
% % %             end     
% % %             
% % %             data = healthy1Data(dataSelect.settings, dataSelect.basePath, dataSelect.patientSelection);
% % %             
% % %         case 'TRI'
% % %             
% % %     end
% % % 
% % %     data.partition;
% % % end
% % % 
% % % function [dimReduct, trainingData, trainingLabel, testingData, testingLabel] = dimReductProcessing(dimReductSelect, data)
% % % 
% % %     switch dimReductSelect.name
% % %         case 'None'
% % %             dimReduct = [];
% % %             
% % %         case 'PCA'
% % %             if ~isfield(dimReductSelect, 'degree')
% % %                 % generate default settings
% % %                 dimReductSelect.degree = 2;
% % %             end
% % %             
% % %             dimReduct = PCATransform(dimReductSelect.degree);
% % %             
% % %         case 'FDA'
% % %             if ~isfield(dimReductSelect, 'degree')
% % %                 % generate default settings
% % %                 dimReductSelect.degree = 2;
% % %             end
% % %             
% % %             dimReduct = FDATransform(dimReductSelect.degree);
% % %     end
% % %     
% % %     if strcmpi(dimReductSelect.name, 'None')
% % %         trainingData = data.trainingData;
% % %         trainingLabel = data.trainingLabel;
% % %         
% % %         testingData = data.testingData;
% % %         testingLabel = data.testingLabel;
% % %     else
% % %         dimReduct.train(data.trainingData, data.trainingLabel);
% % %         
% % %         trainingData = dimReduct.apply(data.trainingData')';
% % %         trainingLabel = data.trainingLabel;
% % %         
% % %         testingData = dimReduct.apply(data.testingData')';
% % %         testingLabel = data.testingLabel;
% % %     end
% % %     
% % % %     figure
% % % %     hold on
% % % %     plot(testingData(1, :), testingData(2, :), '*k')
% % % %     find1 = find(testingLabel == 1);
% % % %     plot(testingData(1, find1), testingData(2, find1), '*r')
% % % %     hold on
% % % %     find0 = find(testingLabel == 0);
% % % %     plot(testingData(1, find0), testingData(2, find0), '*b')
% % % %     grid on
% % % 
% % % end
% % % 
% % % function [classifier, classifierLabel] = classifierProcessing(classifierSelect, trainingData, trainingLabel, testingData)
% % %     switch classifierSelect.name 
% % %         case 'kNN'
% % %             if ~isfield(classifierSelect, 'k')
% % %                 classifierSelect.k = 5;
% % %                 classifierSelect.distanceMetric = 'Euclidean';
% % %             end
% % %             
% % %             classifier = kNNClassifier(classifierSelect.k, classifierSelect.distanceMetric);
% % %             
% % %         case 'QDA'
% % %             
% % %         case 'RBF'
% % %             if ~isfield(classifierSelect, 'nodes')
% % %                 % generate default settings
% % %                 classifierSelect.nodes = 20;
% % %             end
% % %             
% % %             classifier = RBFClassifier(classifierSelect.nodes);
% % %             
% % %         case 'SVM'
% % %             
% % %         case 'NN'
% % %     end
% % %     
% % %     classifier.train(trainingData', trainingLabel);
% % %     classifierLabel = classifier.classify(testingData');
% % % end
% % % 
% % % function [correctTotal, correct, falsePos, falseNeg] = accuracyAssess(classifierLabel, testingLabel)
% % %     correctTotal = sum(testingLabel);
% % %     
% % %     % if they are in both arrays...
% % %     correct = sum(and(classifierLabel, testingLabel));
% % %     
% % %     % false positive
% % %     falsePos = sum(and(xor(classifierLabel, testingLabel), classifierLabel));
% % %     
% % %     % false negative
% % %     falseNeg = sum(and(xor(classifierLabel, testingLabel), testingLabel));
% % % end