classdef SVMIncClassifier < AClassifier
   %This is the Support Vector Machine classifier
   %Uses matlab built in SVM functionality
   properties(Access=public)
       SVMmodel_sup = [];
       SVMmodel_inc = [];
       modelParam = [];
       
       drRetrainingMasterSwitch = 0; % do we want to retrain everything in case DR changes? 
       
       fsmImu = [];
       fsmEkf = [];
       fsmPp = [];
       
       dimReduct = [];
       
       fgModel = [];
       fgStateMachine = [];
       
       gndTruthStateMachine = [];
       
       clusterNormal = [];
       clusterDR = [];
       
       pickPlaceDR = [];
       
       linearSC = [];
       
       % ground truth selector: if sum is more than 1, then voting is used.
       % in ties, the point is discarded
       incsvmFlag = 1;
       incsvmMode = {'highProb'}; % 'highProb_30', 'fsmImu_fsmLabels', 'fsmImu_svmfsmLabels', 'fsmImu_svmLabels'
       incsvmParam = {'30'};
       incsvmUpdate = {0};
       incsvmRatios = {0};
       
       searchArray = 10.^(-2:2:2);
%        searchArray = [];

       % how often do we reset the classifier's training data
       resetTrainingDataAfterClassification = 'exercise'; % never, exercise, patient
       testingDataParam = [];
       
       retrainingRate = 'samples'; % mismatch samples
       retrainingFreq = 50; % '50', 'mismatch10' % when we get this much new 'good' data points, retraing the SVM
       retrainingWindow = 50;
       
       lengthTrainingDataThreshold = 10^6; % if not resetting, this limits the amount of training data entries possible
       lagTimer = 1; % don't incorporate the starting datapoints
       pruneTrainingPoints = 1; % after training, discard all training points that is not part of the SVs
       incLength = 0; % how long to run the inc session for. 0 means run as much as possible. otherwise, it's by sample count
       
       % prob mode: none platt dist1
       svmprobMode = 'platt';

       applyZeroMean = 0; % 1 = remove mean and apply unity covar before training and testing
       downsampleUpperBound = 10000; % the amount of points used for parameter searching
       meanOffset = [];
       stdDevOffset = [];
       
       trainingInputOriginal = {}; % for pick/place and cluster routines
       trainingOutputOriginal = {};
       
       trainingInputOriginalSV = {}; % SV pared down
       trainingOutputOriginalSV = {};
       
       trainingInput = {}; % first cell is the original. subsequent cells are add-ons
       trainingOutput = {};
       
       %Default Kernel is RBF
       kernel = 2;
       degreeFlag = 0;
       gammaFlag = 0;
       coef0Flag = 0;
       
       ignoreVal = -Inf;
    end
    
    methods
        function init(obj,classifierSelect)
            % set up kernel and tuning variables for the SVM training
            
            switch classifierSelect.kernel
                case 'linear'
                    obj.kernel = 0;
                    obj.degreeFlag = 0;
                    obj.gammaFlag = 0;
                    obj.coef0Flag = 0;
                    
                case 'polynomial'
                    obj.kernel = 1;
                    obj.degreeFlag = 1;
                    obj.gammaFlag = 1;
                    obj.coef0Flag = 1;
                    
                case 'radial'
                    obj.kernel = 2;
                    obj.degreeFlag = 0;
                    obj.gammaFlag = 1;
                    obj.coef0Flag = 0;
                    
                case 'sigmoid'
                    obj.kernel = 3;
                    obj.degreeFlag = 0;
                    obj.gammaFlag = 1;
                    obj.coef0Flag = 1;
                    
                otherwise
                    disp('Incorrect input to SVM.init()');
                    disp('Expected Kernel Function ''linear'',''polynomial'',''radial'',''sigmoid''');
            end
                        
            switch classifierSelect.tuning
                case '1'
                    obj.searchArray = 10.^(-2:2:2);
                case '0'
                    obj.searchArray = [];
            end
                               
            switch classifierSelect.resetRate
                case '0'
                    obj.resetTrainingDataAfterClassification = 'never';
                case '1'
                    obj.resetTrainingDataAfterClassification = 'exercise';
                case '2'
                    obj.resetTrainingDataAfterClassification = 'patient';
            end
                                
            switch classifierSelect.retrainRate
                case '1'
                    obj.retrainingRate = 'samples';
                    obj.retrainingFreq = 50;
                case '2'
                    obj.retrainingRate = 'mismatch';
                    obj.retrainingFreq = 0.9;
                case '3'
                    obj.retrainingRate = 'boosting';
                    obj.retrainingFreq = 50;
            end

            obj.incLength = str2num(classifierSelect.incLength);                    
            
            obj.incsvmMode{1} = classifierSelect.adaptMethod;
            param = [];
            switch classifierSelect.adaptMethod
                case 'highProb'
                    param.probThreshold = str2num(classifierSelect.adaptParam{1})/100;
                    
                case {'fsmImu', 'fsmEkf', 'fsmPp'}
                    param.labelSource = classifierSelect.adaptParam{1};
                    
                case 'fgModel'
                    
                case 'clusterNormal'
                    
                case {'clusterDR'}
                    param.dataFormatName = classifierSelect.adaptParam{1};
                    param.dataFormatDofStr = classifierSelect.adaptParam{2};
                    param.dataFormatfitting = classifierSelect.adaptParam{3};
                    
                    param.clusterCountPerClass = str2num(classifierSelect.adaptParam{4}); % number of clusters per class
                    param.dataFormatQuad = str2num(classifierSelect.adaptParam{5});
                    
                    param.updateMachine = str2num(classifierSelect.adaptParam{6}); % update the machine/classifier?
                    
                    param.dataFormatDistanceMetric = classifierSelect.adaptParam{7};
                    param.plateSelection = classifierSelect.adaptParam{8};
                    
                    param.removeMean = classifierSelect.adaptParam{9}; % 'none' is no, 'mean' is remove mean, 'circle' is remove centroid fit
                    param.removeVariance = classifierSelect.adaptParam{10}; % 'none' is no, 'max' 'maxmax' is remove mean, 'circle' is remove centroid fit
                    
                    % percentage of the original data to keep
                    param.percentageOriginalDataToKeepFirstReplacement = str2num(classifierSelect.adaptParam{11})/100; 
                    
                    % percentage of the original data to keep after the first replacement
                    param.percentageOriginalDataToKeepSubsequentReplacement = str2num(classifierSelect.adaptParam{12})/100;
                    
                     % distance from the core to keep in percentage
                    param.percentageDistanceFromCoreToKeep = str2num(classifierSelect.adaptParam{13})/100;

                case 'pickPlaceDR'
                    param.datapointsInACluster = str2num(classifierSelect.adaptParam{1}); % number of data points for each cluster
                    param.updateMachine = str2num(classifierSelect.adaptParam{2}); % update the machine/classifier?
                    
                    % overlap percentage allowed in distance
                    param.overlapInDistance = str2num(classifierSelect.adaptParam{3})/100;

                    % age of old datapoints
                    param.ageThresholdForOldDatapoints = str2num(classifierSelect.adaptParam{4});
                    
                case 'groundTruth'
                    param.labelSource = classifierSelect.adaptParam{1};
                    
                otherwise
                    % use default settings
            end
            
            obj.incsvmParam{1} = param;
        end
        
        function initializeAdaptiveElements(obj)
            segmentInfo.startTimeVal =[];
            segmentInfo.startTimeInd = [];
            segmentInfo.endTimeVal =[];
            segmentInfo.endTimeInd = [];
            segmentInfo.segmentName =[];
            segmentInfo.segmentInclude = [];
            
            obj.fgStateMachine.startInd = 1;
            obj.fgStateMachine.previousSegmentInfo = segmentInfo;
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
            
            if ~isempty(obj.pickPlaceDR)
                obj.pickPlaceDR.updateDimReductClass(dimReduct);
            end
        end
        
        function error = train(obj,input,output)
            % this function performs a simple parameter search, then trains
            % the SVM based on the optimal parameter
            
            % save the original training data for cluster and pick/place
            % side classifiers, since the SVs may not be the best
            % datapoints to train these side classifiers
            obj.trainingInputOriginal = input;
            obj.trainingOutputOriginal = output;
            
            % update the internal library for inc svm 
            obj.trainingInput{1} = input;
            obj.trainingOutput{1} = output;
            
            % train the initial classifier with the supervised data
            obj.SVMmodel_sup = incrementalTraining(obj);
            
            % save the original training SVs on the side for when we need
            % to reset the classifiers, after incrementalTraining have
            % removed non-SVs from the dataset
            obj.trainingInputOriginalSV = obj.trainingInput;
            obj.trainingOutputOriginalSV = obj.trainingOutput;
            
% %             obj.plotSupportVectors(inputUnity, output, obj.SVMmodel);
            
            inputDR = preprocessDataTraining(obj, input);

            [label_training, error, decisionValue] =...
                svmpredict_mex(output, inputDR, obj.SVMmodel_sup, '-q');
            error = 1-error(1)/100;
            
              if 0
                  figure;
                  plot(decisionValue, '.');
                  hold on
                  svIndices = obj.SVMmodel_sup.sv_indices;
                  plot(svIndices, decisionValue(svIndices), 'ro');
                  plot(label_training, 'x');
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
            % reset the training sets. the classify function should call
            % the retrain command after this though
            obj.trainingInput = obj.trainingInputOriginalSV;
            obj.trainingOutput = obj.trainingOutputOriginalSV;
%             obj.SVMmodel_inc = incrementalTraining(obj);
            
            % reset the side classifiers. the other side classifiers are 
            % static
            % then update the side classifiers too
            if ~isempty(obj.clusterDR)
                obj.clusterDR.resetTrainingData();
            end
            
            if ~isempty(obj.pickPlaceDR)
                obj.pickPlaceDR.resetTrainingData();
            end
        end
               
        function [out_inc, prob_inc] = classify(obj,input,extra,incMode)
            % NOTE THAT 'INPUT' INTO THIS FUNCTION IS NOT DR'ed, until
            % after the first line, so there's two copies of this data now
            retrainFlag = 0;
            
            inputDR = preprocessDataTesting(obj, input);
                        
            if obj.incsvmFlag && incMode % set at top, and in the function passed in
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
                
                % retrain the SVM. the dim reduct may have changed, or the
                % training data contents may have all be rearranged due to
                % resetting
                fprintf('New exercise SVM refresh training ');
                obj.SVMmodel_inc = incrementalTraining(obj);
                fprintf('\n');
                
                if obj.drRetrainingMasterSwitch
                    % if the side classifiers uses cluster/pickplace, we should
                    % reset those based on any changes of DR as well, keeping
                    % all the existing data
                    if ~isempty(obj.clusterDR)
                        obj.clusterDR.retrain([], 1, -1);
                    end
                    
                    if ~isempty(obj.pickPlaceDR)
                        obj.pickPlaceDR.retrain([],[]);
                    end
                end
                
                out_inc = -1*ones(size(inputDR, 1), 1);
                prob_inc = -1*ones(size(inputDR, 1), 1);
                decisionValue_inc = -1*ones(size(inputDR, 1), 1); 
                
                newInput0 = [];  newInput1 = [];  newDV0 = [];
                newOutput0 = []; newOutput1 = []; newDV1 = [];
                updateTime = []; 
                
                % set an ending time for the adaptiveness
                if obj.incLength == 0
                    incEndTime = size(inputDR, 1)-obj.retrainingWindow;
                else
                    % select the time corresponding to the end of a given
                    % segment
                    currSegSelect = obj.incLength;
                    if currSegSelect > length(extra.useTime.endTimeVal)
                        currSegSelect = length(extra.useTime.endTimeVal);
                    end
                    
                    % pull out the timestamp that corresponds with the end
                    % of the segment of interest, plus 0.25 seconds steps into the
                    % future
                    [~, incEndTime] = findClosestValue(extra.useTime.endTimeVal(currSegSelect) + 0.25, extra.testingTime);
                end
                
                if incEndTime > size(inputDR, 1)
                    incEndTime = size(inputDR, 1)-obj.retrainingWindow;
                end

                for ind_input = 1:obj.retrainingWindow:incEndTime
                    % classify
                    ind_range = ind_input:ind_input+obj.retrainingWindow;
                    currentInput = input(ind_range, :);
                    currentInputDR = inputDR(ind_range, :);
                    
                    [out_inc(ind_range), ~, decisionValue_inc(ind_range)] = ...
                        svmpredict_mex(zeros(size(currentInputDR,1),1), currentInputDR, obj.SVMmodel_inc, '-q');
                    
%                     [out_inc2, ~, decisionValue_inc2] = svmpredict_mex(zeros(size(inputDR,1),1), inputDR, obj.SVMmodel_inc, '-q');

                    % select which data points to incorporate into inc svm
                    [ind0loc, ind1loc, out_sideClass] = ...
                            incSvmIndSelection(obj, out_inc, decisionValue_inc, currentInput, extra, ind_range, retrainFlag);

%                     if length(ind0locArray) == 1 
%                         % only one entry in the array
%                        ind0loc = ind0locArray{1};
%                        ind1loc = ind1locArray{1};
%                     else
%                         % setxor
%                         ind0loc = ind_range;
%                         ind1loc = ind_range;
%                         for ind_setxor = 1:length(ind0locArray)
%                             ind0loc = setxor(ind0loc, ind0locArray{ind_incsvmmode});
%                             ind1loc = setxor(ind1loc, ind1locArray{ind_incsvmmode});
%                         end
%                     end
                    
                    % accumulated enough points. add to the training, and
                    % retrain, or retrain after fgmodel declares a new
                    % segment
                    retrainFlag = 0;
                    
                    % adding the samples-to-retrain with to our buffer
                    switch obj.retrainingRate
                        case {'samples', 'mismatch'}
                            % add in the new obs candidates for svm
                            % training. fill the retraining array as they
                            % come
                            if ind_range(end) > obj.lagTimer
                                if ~isempty(ind0loc) %&& length(newOutput0) < obj.retrainingFreq
                                    newInput0 = [newInput0; input(ind0loc, :)];
                                    newOutput0 = [newOutput0; zeros(size(ind0loc))'];
                                end
                                
                                if ~isempty(ind1loc) %&& length(newOutput1) < obj.retrainingFreq
                                    newInput1 = [newInput1; input(ind1loc, :)];
                                    newOutput1 = [newOutput1; ones(size(ind1loc))'];
                                end
                            end
                            
                        case 'boosting'
                            if ind_range(end) > obj.lagTimer
                                % add in the new obs candidates for svm
                                % training. fill the array only if they are
                                % close to the margin (ie will change the
                                % margins of the SV distro)
                                if ~isempty(ind0loc)
                                    % test the ind0loc points. these points
                                    % should be 0 (according to the side
                                    % classifier) but is not
                                    samplePointsDR = inputDR(ind0loc, :);
                                    [out_boostTest, ~, decisionValue_boostTest] = ...
                                        svmpredict_mex(zeros(size(samplePointsDR,1),1), samplePointsDR, obj.SVMmodel_inc, '-q');
                                    
                                    decisionValuePass = decisionValue_boostTest > -1;
                                    ind0locKeep = ind0loc(decisionValuePass);
                                    newInput0 = [newInput0; input(ind0locKeep, :)];
                                    newOutput0 = [newOutput0; zeros(size(ind0locKeep))'];
                                    newDV0 = [newDV0; decisionValue_boostTest(decisionValuePass)];
                                end
                                
                                if ~isempty(ind1loc)
                                    samplePointsDR = inputDR(ind1loc, :);
                                    [out_boostTest, ~, decisionValue_boostTest] = ...
                                        svmpredict_mex(zeros(size(samplePointsDR,1),1), samplePointsDR, obj.SVMmodel_inc, '-q');
                                    
                                    decisionValuePass = decisionValue_boostTest < 1;
                                    ind1locKeep = ind1loc(decisionValuePass);
                                    newInput1 = [newInput1; input(ind1locKeep, :)];
                                    newOutput1 = [newOutput1; ones(size(ind1locKeep))'];
                                    newDV1 = [newDV1; decisionValue_boostTest(decisionValuePass)];
                                end
                            end
                    end
                    
                    % is the buffer sufficiently full? we can update the
                    % svm if it is
                    switch obj.retrainingRate
                        case {'samples', 'boosting'}
                            if (length(newOutput0) >= obj.retrainingFreq && ...
                                    length(newOutput1) >= obj.retrainingFreq) ...
                                    || ...
                                    (~isempty(newOutput0) && ~isempty(newOutput1) && ...
                                    length(obj.incsvmMode) == 1 && strcmp(obj.incsvmMode{1}, 'fgModel'))
                                % retrain either when the buffer is full,
                                % or a segment is detected for fgmodel and
                                % fg model is the only side classifier
                                % available
                                retrainFlag = 1;
                                
                                if length(newOutput0) >= obj.retrainingFreq && length(newOutput1) >= obj.retrainingFreq
                                    lengthToUse = obj.retrainingFreq;
                                else
                                    lengthToUse = min([length(newOutput0) length(newOutput1)]);
                                end
                            end
                            
                        case 'mismatch'                            
                            if sum(out_sideClass == out_inc(ind_range)) / obj.retrainingWindow < obj.retrainingFreq && ...
                                    (length(newOutput0) >= obj.retrainingWindow && ...
                                    length(newOutput1) >= obj.retrainingWindow) ...
                                    
                                retrainFlag = 1;
                                
                                if length(newOutput0) >= obj.retrainingWindow && length(newOutput1) >= obj.retrainingWindow
                                    lengthToUse = obj.retrainingWindow;
                                else
                                    lengthToUse = min([length(newOutput0) length(newOutput1)]);
                                end
                            end
                    end
                    
                    if retrainFlag
                        % add new block to the SVM training set
                        if ~isempty(newDV1)
                            % pull out the points that are furthest from
                            % the margin
                            [blah0, newInput0IndToUseTemp] = sort(newDV0, 'descend');
                            [blah1, newInput1IndToUseTemp] = sort(newDV1, 'ascend');
                            
                            newInput0IndToUse = newInput0IndToUseTemp(1:lengthToUse);
                            newInput1IndToUse = newInput1IndToUseTemp(1:lengthToUse);
                        else
                            % just use the first so and so
                            newInput0IndToUse = 1:lengthToUse;
                            newInput1IndToUse = 1:lengthToUse;
                        end
                        
                        newInput = [newInput0(newInput0IndToUse, :); newInput1(newInput1IndToUse, :)];
                        newOutput = [newOutput0(newInput0IndToUse, :); newOutput1(newInput1IndToUse, :)];
                        
                        newTrainingIndex = length(obj.trainingOutput) + 1;
                        obj.trainingInput{newTrainingIndex} = newInput;
                        obj.trainingOutput{newTrainingIndex} = newOutput;
                        
                        % clear out the cache
                        newInput0 = [];  newInput1 = [];  newDV0 = [];
                        newOutput0 = []; newOutput1 = []; newDV1 = [];
                        
                        % retrain the SVM
                        fprintf('Updating at %u/%u (inc cap at %u)\n', ind_input, size(inputDR, 1), incEndTime);
                        obj.SVMmodel_inc = incrementalTraining(obj);
                        fprintf('\n');
                        
                        updateTime = [updateTime ind_input];
                    end
                end
                
                % there may be obs entries cut off by the inc testing
                % array. run the rest normally
                startInd = find(out_inc == -1, 1, 'first');
                
                if ~isempty(startInd)
                    % classify
                    ind_range = startInd:size(inputDR, 1);
                    currentInputDR = inputDR(ind_range, :);
                    [out_inc(ind_range), ~, decisionValue_inc(ind_range)] = ...
                        svmpredict_mex(zeros(size(currentInputDR,1),1), currentInputDR, obj.SVMmodel_inc, '-q');
                    prob_inc(ind_range) = sigmoid_predict(obj, decisionValue_inc(ind_range), obj.SVMmodel_inc.ProbA, obj.SVMmodel_inc.ProbB);
                end

            else
                % no adaptiveness. this likely is reached by the training
                % component, and not by the testing 
                [out_sup, err_sup, decisionValue_sup] = ...
                    svmpredict_mex(extra.trainingLabel, inputDR, obj.SVMmodel_sup, '-q');                
                
                prob_sup = sigmoid_predict(obj, decisionValue_sup, obj.SVMmodel_sup.ProbA, obj.SVMmodel_sup.ProbB);
                
                out_inc = out_sup;
                prob_inc = prob_sup;
                updateTime = [];
            end
            
            if 0
                if exist('extra', 'var')
                    % for the original data input
                    graphTime = extra.time;
                    graphData = extra.data(:, 1:5);
                    segStartTime = extra.testingSegTime.startTimeVal;
                    segEndTime = extra.testingSegTime.endTimeVal;
                    
                    % for the classifier time
                    testingTime = extra.testingTime;
                    titleStr = extra.identString;
                    
                    updateTime = testingTime(updateTime);
                else
                    graphTime = 1:length(prob_sup);
                    graphData = input(:, 1:5);
                    
                    testingTime = 1:length(prob_sup);
                    segStartTime = [];
                    segEndTime = [];
                    titleStr = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
                end
                
                % plot the stuff
                h =    figure;
                axes(1) = subplot(311);
                plot(graphTime, graphData);
                plotBoxes(h, segStartTime, segEndTime, 'b', 0);   
                title(titleStr);
  
                axes(2) = subplot(312);
                localAcc = sum(out_sup == labels)/length(labels);
                hold on
                plot(testingTime, prob_sup, 'r');
                plot(testingTime, out_sup, 'ro'); 
                plot([testingTime(1) testingTime(end)], [0.5 0.5], 'k');
                plot([testingTime(1) testingTime(end)], [0.5 0.5]+0.3, 'k');
                plot([testingTime(1) testingTime(end)], [0.5 0.5]-0.3, 'k');
                title(['sup svm - classifier (r) and ground truth (b). Acc = ' num2str(localAcc)]);
                if exist('labels', 'var')
                    plot(testingTime, labels, 'bx');
                end
                ylim([-0.1 1.1]);
                
                axes(3) = subplot(313);
                localAcc = sum(out_inc == labels)/length(labels);
                hold on
                plot(testingTime, prob_inc, 'r');
                plot(testingTime, out_inc, 'ro'); 
                plot([testingTime(1) testingTime(end)], [0.5 0.5], 'k');
                plot([testingTime(1) testingTime(end)], [0.5 0.5]+0.3, 'k');
                plot([testingTime(1) testingTime(end)], [0.5 0.5]-0.3, 'k');
                title(['inc svm ' num2str(incMode) ' - classifier (r) and ground truth (b). Acc = ' num2str(localAcc)]);
                if exist('labels', 'var')
                    plot(testingTime, labels, 'bx');
                end
                plotBoxes(h, updateTime, updateTime, 'g', 0); 
                ylim([-0.1 1.1]);
                
                linkaxes(axes, 'x');
                
                saveas(h, ['D:\MATLABResults\Seg_Healthy1\figs\IncSVM-' titleStr '-' num2str(localAcc) '-' datestr(now, 'HH-MM-SS') '.fig']);
                close(h);
            end
        end
        
        function adaptiveComponent(obj, trainingData)
            % train the various different adaptive components that would be
            % used 
            
            for i = 1:length(obj.incsvmMode)
                switch obj.incsvmMode{i}
                    case 'highProb'
                        
                    case 'fsmImu'
                        auxStateImu = trainingData.incTraining.auxStateImu;
                        obj.fsmImu = finiteStateMachine();
                        obj.fsmImu.train(auxStateImu);
                        
                    case 'fsmEkf'
                        auxStateEkf = trainingData.incTraining.auxStateEkf;
                        obj.fsmEkf = finiteStateMachine();
                        obj.fsmEkf.train(auxStateEkf);
                        
                    case 'fsmPp'
                        auxStatePp = trainingData.incTraining.auxStatePp;
                        obj.fsmPp = finiteStateMachine();
                        obj.fsmPp.train(auxStatePp);
                        
                    case 'fgModel'
                        obj.fgModel = trainingData.incTraining.fgModel;
                        
                    case 'clusterNormal'
                        
                    case 'groundTruth'
                        
                    case 'clusterDR'
                        obj.incsvmParam{i}.dataFormatDofNum = trainingData.sigDofArray;
                        obj.trainClusterDR(obj.incsvmParam{i}, trainingData);
                        
                    case 'pickPlaceDR'
                        obj.trainPickPlaceDR(obj.incsvmParam{i});
                end
            end
        end
        
        function trainClusterDR(obj, incsvmParam, trainingData)                
            trainingDataParse = clusterDataGatewayTraining(trainingData, incsvmParam);
            
            obj.clusterDR = clusterMachine(incsvmParam);
            obj.clusterDR.updateDimReductClass(obj.dimReduct); % apply this now because this classifier would not have existed yet when the DR update for the SVM is happening
            obj.clusterDR.retrain(trainingDataParse, 1, incsvmParam.percentageDistanceFromCoreToKeep); % DR is performed in the training function
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