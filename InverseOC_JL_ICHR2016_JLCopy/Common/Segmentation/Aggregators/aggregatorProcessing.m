function [classifier, classifierLabel] = ...
    aggregatorProcessing(aggregatorSelect, baseClassifier, baseClassifierLabel, trainingData, trainingLabel)
    % perform the aggregator and classification training
    
    switch aggregatorSelect.name
        case 'None'
            % no aggregator. pass the parameters back out
            classifier = baseClassifier;
            classifierLabel = baseClassifierLabel;
            
        case 'Boosting'
            % combine 3 base algorithms in an boosting algorithm
            classifier = BoostAggregator(baseClassifier);
            if isfield(aggregatorSelect, 'iteration')
                classifier.init(aggregatorSelect.iteration);
            end

            classifier.train(trainingData, trainingLabel);
            classifierLabel = classifier.classify(trainingData); 
            
        case 'Bagging'
             classifier = BaggingAggregator(baseClassifier);
            if isfield(aggregatorSelect, 'iteration')
                classifier.init(aggregatorSelect.iteration);
            end
             
            classifier.train(trainingData, trainingLabel);
            classifierLabel = classifier.classify(trainingData); 
    end
end