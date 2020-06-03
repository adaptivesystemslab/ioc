classdef clusterClassifier < AClassifier
   %A wrapper to conver the cluster machine into the cluster classifier
   properties(Access=public)
       clusterDR = [];
       
       incsvmParam = [];
       dimReduct = [];
       
       retrainingWindow = 50;
       resetTrainingDataAfterClassification = 'exercise';
    end
    
    methods
        function obj = clusterClassifier(classifierSelect)
            % constructor
            obj.init(classifierSelect);
        end
        
        function init(obj, classifierSelect)
            
            param.dataFormatName = classifierSelect.adaptParam{1}(1:2);
            param.dataFormatDof = str2num(classifierSelect.adaptParam{1}(3));
            param.clusterCountPerClass = str2num(classifierSelect.adaptParam{2}); % number of clusters per class
            param.dataFormatQuad = str2num(classifierSelect.adaptParam{3});
            
            param.updateMachine = str2num(classifierSelect.adaptParam{4}); % update the machine/classifier?
            
            % percentage of the original data to keep
            param.percentageOriginalDataToKeepFirstReplacement = str2num(classifierSelect.adaptParam{5})/100;
            
            % percentage of the original data to keep after the first replacement
            param.percentageOriginalDataToKeepSubsequentReplacement = str2num(classifierSelect.adaptParam{6})/100;
            
            % distance from the core to keep in percentage
            param.percentageDistanceFromCoreToKeep = str2num(classifierSelect.adaptParam{7})/100;
            
            obj.incsvmParam = param;
        end
        
        function out = preprocessDataTraining(obj, input)
            % apply the dim reduction 
            out =  obj.dimReduct.apply(input);
            
            % for SVM, the data should be normalized before inserting into
            % the system: zero mean and unity covariance. store the offset
            % covar to apply to the testing data
            if obj.applyZeroMean
                [out, obj.meanOffset, obj.stdDevOffset] = zeroMeanCovarSignal(out, 1, 1);
            end
        end
        
        function out = preprocessDataTesting(obj, input)
            % apply the dim reduction
            out =  obj.dimReduct.apply(input);
            
            if obj.applyZeroMean
                % shift the data to 0 mean, 1 covar, if applicable, from
                % previously generated offset vectors
                out = zeroMeanCovarSignal(out, 1, 1, obj.meanOffset, obj.stdDevOffset);
            end
        end
        
        function updateDimReductClass(obj, dimReduct)
            % save the dim reduction into the class
            obj.dimReduct = dimReduct;
            
            % then update the side classifiers too
            if ~isempty(obj.clusterDR)
                obj.clusterDR.updateDimReductClass(obj.dimReduct);
            end
        end
        
        function [input, output] = combineDataForTraining(obj)
            input = [];
            output = [];
            
            % combine all the collected data so far
            for ind = 1:length(obj.trainingInput)
                input = [input; obj.trainingInput{ind}];
                output = [output; obj.trainingOutput{ind}];
            end
        end
        
        function svmModel = incrementalTraining(obj)
            % concat all the training data
            [input, output] = combineDataForTraining(obj);
            
            % apply the dim reduction
            inputSVMTrain = preprocessDataTraining(obj, input);
            
            if 0
                p1points = inputSVMTrain(output == 1, :);
                p1ds = 1:size(p1points, 1);
                
                p0points = inputSVMTrain(output == 0, :);
                p0ds = 1:size(p0points, 1);
                
                scatter3(p1points(p1ds, 1), p1points(p1ds, 2), p1points(p1ds, 3), 'bo')
                hold on; title('')
                scatter3(p0points(p0ds, 1), p0points(p0ds, 2), p0points(p0ds, 3), 'rx')
            end
            
            fprintf(' (SVM update w/ %u dp) ', length(output));
            
            % if applicable, tune for SVM parameters
            newModelParam = SVMClassifier_svmTuning(obj, inputSVMTrain, output, obj.searchArray);
            
            % now for the SVM augmentations
            switch obj.svmprobMode
                case 'none'
                case 'platt'
                    newModelParam.string = [newModelParam.string ' -b 1'];
                case 'dist1'
            end
            
            % perform the actual testing
            svmModel = svmtrain_mex(output, inputSVMTrain, newModelParam.string);
            
            if obj.pruneTrainingPoints
                % keep only the SVs (in original data form, not dim reduct)
                svIndices = svmModel.sv_indices;
                svInput = input(svIndices, :);
                svOutput = output(svIndices);

                % erase the current training set
                obj.trainingInput = {};
                obj.trainingOutput = {};
                
                % and save the ones we want to keep
                obj.trainingInput{1} = svInput;
                obj.trainingOutput{1} = svOutput;
            end
        end
        
        function resetTheClassifier(obj)
            % reset the side classifiers. the other side classifiers are 
            % static
            % then update the side classifiers too
            if ~isempty(obj.clusterDR)
                obj.clusterDR.resetTrainingData();
            end
        end
        
        function inputCluster = getInput(obj, extra, ind_range)
            switch obj.clusterDR.dataFormatName
                case 'PC'
                    if ~isempty(extra.testingData)
                        inputCluster = extra.testingData(ind_range, :);
                    else
                        inputCluster = extra.trainingData(ind_range, :);
                    end
                    
                case 'PP'
                    inputCluster = extra.data(ind_range, :);
            end
        end
               
        function [out_class, prob_class] = classify(obj,input,extra,incMode)
            
            switch obj.resetTrainingDataAfterClassification
                case 'never'
                    % never reset the classifier
                    
                case 'exercise'
                    % start each exercise fresh
                    fprintf('Resetting incSVM (exercise)... ');
                    resetTheClassifier(obj);
                    
                case 'patient'
                    currDataParam.subject = extra.subjectData{1}.subjectNumber;
                    currDataParam.session = extra.subjectData{1}.sessionNumber;
                    currDataParam.exercise = extra.subjectData{1}.exerciseName;
                    
                    if isempty(obj.testingDataParam)
                        % first run. don't need to reset
                        
                    elseif obj.testingDataParam.subject == extra.subjectData{1}.subjectNumber
                        % update the param but don't reset
                        
                    else
                        % they differ. reset the classifier
                        fprintf('Resetting incSVM (patient)... ');
                        resetTheClassifier(obj);
                    end
                    
                    % update the subject data that we're looking at
                    obj.testingDataParam = currDataParam;
                    
                otherwise
                    % we probably should prune it anyway. if there is more than
                    % 'lengthTrainingDataThreshold' sets of data...
                    if length(obj.trainingInput) > obj.lengthTrainingDataThreshold
                        valuesToKeep = randperm(length(obj.trainingOutput), obj.lengthTrainingDataThreshold);
                        resetInput = obj.trainingInput(valuesToKeep, :);
                        resetOutput = obj.trainingOutput(valuesToKeep, :);
                        obj.trainingInput = resetInput;
                        obj.trainingOutput = resetOutput;
                        obj.SVMmodel_inc = incrementalTraining(obj);
                    end
            end
            
            
            
            out_class = -1*ones(size(input, 1), 1);
            prob_class = -1*ones(size(input, 1), 1);
            
            incEndTime = size(input, 1)-obj.retrainingWindow;
            
            for ind_input = 1:obj.retrainingWindow:incEndTime
                ind_range = ind_input:ind_input+obj.retrainingWindow;
                inputCluster = getInput(obj, extra, ind_range);
                
                [out_class(ind_range), prob_class(ind_range)] = obj.clusterDR.classify(inputCluster);
                
                if obj.incsvmParam.updateMachine && incMode
                    % how many of the old data points to keep
                    percentageToKeep = obj.incsvmParam.percentageOriginalDataToKeepSubsequentReplacement;
                    
                    trainingDataParse.input = inputCluster; % pre-PCA data
                    trainingDataParse.output = out_class(ind_range);
                    trainingDataParse.distance = prob_class(ind_range);
                    
                    if ~isempty(extra.testingData)
                        testingLabel = extra.testingLabel(ind_range);
                    else
                        testingLabel = extra.trainingLabel(ind_range);
                    end
                    
                    trainingDataParse.outputGndTruth = testingLabel; % INDEXING
                    
                    obj.clusterDR.retrain(trainingDataParse, percentageToKeep, obj.incsvmParam.percentageDistanceFromCoreToKeep);
                end
                
            end
            
            startInd = find(out_class == -1, 1, 'first');
            
            if ~isempty(startInd)
                % classify
                ind_range = startInd:size(input, 1);
                inputCluster = getInput(obj, extra, ind_range);
                
                [out_class(ind_range), prob_class(ind_range)] = obj.clusterDR.classify(inputCluster);
            end
        end
        
        function train(obj, trainingData)
            % train the various different adaptive components that would be
            % used 
            
             obj.trainClusterDR(obj.incsvmParam, trainingData);
        end
        
        function trainClusterDR(obj, incsvmParam, trainingData)
            switch incsvmParam.dataFormatName
                case 'PC'
                    trainingDataParse.input = trainingData.trainingData; % pre-PCA data
                    trainingDataParse.output = trainingData.trainingLabel;
                    trainingDataParse.outputGndTruth = trainingData.trainingLabel;
                case 'PP'
                    trainingDataParse.input = trainingData.data; % the original data
                    trainingDataParse.output = trainingData.label;
                    trainingDataParse.outputGndTruth = trainingData.label;
            end
            
            obj.clusterDR = clusterMachine(incsvmParam);
            obj.clusterDR.updateDimReductClass(obj.dimReduct); % apply this now because this classifier would not have existed yet when the DR update for the SVM is happening
            obj.clusterDR.train(trainingDataParse); % DR is performed in the training function
        end
        
        function trainPickPlaceDR(obj, incsvmParam)
%             sideClassInput = obj.dimReduct.apply(obj.trainingInputOriginal);
%             sideClassOutput = obj.trainingOutputOriginal;
            
            sideClassInput = obj.trainingInputOriginal;
            sideClassOutput = obj.trainingOutputOriginal;

            obj.pickPlaceDR = pickPlaceMachine(incsvmParam);
            obj.pickPlaceDR.updateDimReductClass(obj.dimReduct);
            obj.pickPlaceDR.train(sideClassInput, sideClassOutput);
        end
    
        
        function emission = sigmoid_predict(obj, decision_value, A, B)
            if ~isempty(A) && ~isempty(B)
                fApB = decision_value*A+B;
            else
                fApB = ones(size(decision_value));
            end
            
%             ind = fApB >= 0;  % 1-p used later; avoid catastrophic cancellation
%             emission(ind) = exp(-fApB(ind))./(1+exp(-fApB(ind)));
%             
%             ind = fApB < 0;
%             emission(ind) = 1/(1+exp(fApB(ind)));
            
%             fullEmission1 = exp(-fApB)./(1+exp(-fApB));
%             fullEmission2 = 1/(1+exp(fApB));

            emission = exp(-fApB)./(1+exp(-fApB));
        end

        function obj_copy = copy(obj)
           obj_copy = SVMClassifier();
           kernels = {'linear','polynomial','radial','sigmoid'};
           obj_copy.init(kernels{obj.kernel+1}); % linear would = 0 if the +1 didn't exist
        end
        
        function counter = labelChecker(obj, counter, currLabel, prevLabel)
            % currLabel = the label that was just assessed
            % prevLabel = the label that was before that
            
            % assume previous label are 1 if at the start of the array
            if currLabel == prevLabel
                % reset the counter when the label switches
                counter = counter + 1;
            else
                counter = 1;
            end
        end
    end

end