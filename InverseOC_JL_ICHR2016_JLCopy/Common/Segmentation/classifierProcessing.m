function [classifier, classifierLabel] = classifierProcessing(classifierSelect, aggregatorSelect, trainingDataStruct, dimReduct)
    switch classifierSelect.name 
        case 'kNN'
            if ~isfield(classifierSelect, 'k')
                classifierSelect.k = 5;
                classifierSelect.distanceMetric = 'Euclidean';
            end
            
            classifier = kNNClassifier(classifierSelect.k, classifierSelect.distanceMetric);
          
        case 'LDA'
            classifier = LDAClassifier;
            
        case 'QDA'
            classifier = QDAClassifier;
            
        case 'RBF'
            if ~isfield(classifierSelect, 'nodes')
                % generate default settings
                classifierSelect.nodes = 20;
            end
            
            classifier = RBFClassifier(classifierSelect.nodes);
%             trainingData = trainingData';
            
        case 'SVM'
            classifier = SVMClassifier;
            if isfield(classifierSelect, 'kernel')
                classifier.init(classifierSelect);
            end
            
        case 'SVMInc'
            classifier = SVMIncClassifier;
            if isfield(classifierSelect, 'kernel')
                classifier.init(classifierSelect);
            end
                        
        case 'SVMOneClass'
            classifier = SVMOneClassClassifier;
            if isfield(classifierSelect, 'kernel')
                classifier.init(classifierSelect.kernel);
            end
            
        case 'SVMThreeClass'
            classifier = SVMThreeClassClassifier;
            if isfield(classifierSelect, 'kernel')
                classifier.init(classifierSelect.kernel);
            end
            
        case 'SVMRetrain'
            classifier = SVMRetrainClassifier;
            if isfield(classifierSelect, 'kernel')
                classifier.init(classifierSelect);
            end
            
        case 'NN'
            if ~isfield(classifierSelect, 'net_params')
                % generate default settings
                classifierSelect.net_params = [10 10 10];
            end
            classifier = NNClassifier(classifierSelect.net_params);
            
        case 'SL'
            classifier = SLClassifier;
            
        case 'Threshold'
            classifier = ThresClassifier(classifierSelect.percentage);

        case 'Compose'
            classifier = ComposeClassifier(classifierSelect);
            
        case 'Cluster'
            classifier = clusterClassifier;
            if isfield(classifierSelect, 'kernel')
                classifier.init(classifierSelect);
            end
            
        case 'PolyPCA'
            classifier = PolyPCAClassifier;
    end
    
    if strcmpi(aggregatorSelect.name, 'None')
        if isa(classifier, 'SVMIncClassifier')
            % for the incremental, all the dim transform is done inside the
            % classifier
            classifier.updateDimReductClass(dimReduct);
            classifier.train(trainingDataStruct.trainingData, trainingDataStruct.trainingLabel);
            classifierLabel = classifier.classify(trainingDataStruct.trainingData, trainingDataStruct, 0);
            
        elseif isa(classifier, 'clusterClassifier')
            % classifier
            classifier.updateDimReductClass(dimReduct);
            classifier.train(trainingDataStruct.trainingData, trainingDataStruct.trainingLabel);
            classifierLabel = [];
            
        elseif isa(classifier, 'PolyPCAClassifier')
            classifier.train(trainingDataStruct);
            classifierLabel = [];

        else
            % if no aggregator is being called, train the classifier
            classifier.train(trainingDataStruct.trainingDataDR, trainingDataStruct.trainingLabel);
            classifierLabel = classifier.classify(trainingDataStruct.trainingDataDR);
        end
        
    else
        % otherwise, don't, to save computational cycles, since the
        % aggregator will perform the training
        classifierLabel = [];
    end
end