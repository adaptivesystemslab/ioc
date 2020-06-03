   classdef PolyPCAClassifier < AClassifier
   %This is the Support Vector Machine classifier
   %Uses matlab built in SVM functionality
   
    properties(Access=public)
       dimReductModels = {};
       classifierModels = {};
       motionSet = {};
       
       dataStr = '15_10_1_1101_110';
       dimReductStr = 'PCA_0';
       classifierStr = 'Cluster_radial_0_1_3_0_clusterDR_CR0_2_4_0_50_-1_-1';
       
       applyZeroMean = 0; % 1 = remove mean and apply unity covar before training and testing
       downsampleUpperBound = 10000; % the amount of points used for parameter searching
       meanOffset = [];
       stdDevOffset = [];
       
       %Default Kernel is RBF
       kernel = 2;
       degreeFlag = 0;
       gammaFlag = 0;
       coef0Flag = 0;
       
       ignoreVal = -Inf;
       
       searchArray = 10.^(-2:2:2);
%        searchArray = [];
    end
    
    methods
        function init(obj)
            % set up kernel and tuning variables for the SVM training
            
          
        end
        
        function adaptiveComponent(obj, trainingData)
            
        end
        
        function initializeAdaptiveElements(obj)
%             segmentInfo.startTimeVal =[];
%             segmentInfo.startTimeInd = [];
%             segmentInfo.endTimeVal =[];
%             segmentInfo.endTimeInd = [];
%             segmentInfo.segmentName =[];
%             segmentInfo.segmentInclude = [];
%             
%             obj.fgStateMachine.startInd = 1;
%             obj.fgStateMachine.previousSegmentInfo = segmentInfo;
        end
        
        function input = preprocessData(obj, input)
%             % exp: convert to cylindrical
%             [THETA,RHO,Z] = cart2pol(input(:, 1), input(:, 2), input(:, 3));
%             input = [THETA, RHO, Z];

            % convert to polar coordinates
        end
        
        function [dataSettingsSelect, dimReductSelect, classifierSelect] = setSettings(obj)
            dataSettingsSelect = pollDataSettings(obj.dataStr); 
            dimReductSelect = pollDimReduction(obj.dimReductStr);
            classifierSelect = pollClassifier(obj.classifierStr);
        end
        
        function [error, labels] = train(obj,trainingStruct,output)
            % this function performs a simple parameter search, then trains
            % the SVM based on the optimal parameter
           [dataSettingsSelect, dimReductSelect, classifierSelect] = setSettings(obj);

            % pull out all the motions that are the same
            for ind_motionSet = 1:length(trainingStruct.motionSet)
                motionSetEntry{ind_motionSet} = trainingStruct.motionSet{ind_motionSet}(1:8);
                
                currMotion = trainingStruct.motionSet{ind_motionSet};
                currTrainingData = [];
                currTrainingLabel = [];
                
                for ind_subjData = 1:length(trainingStruct.subjectData)
                    currSubj = trainingStruct.subjectData{ind_subjData};
                    
                    if strcmpi(currMotion(1:8), currSubj.exerciseType(1:8))
                        % is the same
                        [testingData, testingLabel] = qdq2stacked(obj, currSubj, dataSettingsSelect, trainingStruct.subjectData(ind_subjData));
                        
                        currTrainingData = [currTrainingData; testingData];
                        currTrainingLabel = [currTrainingLabel; testingLabel];
                    else
                        continue
                    end
                end
                
                % train the dim reduct algorithm
                currTrainingDataStruct.trainingData = currTrainingData;
                currTrainingDataStruct.trainingLabel = currTrainingLabel;
                obj.dimReductModels{ind_motionSet} = ...
                    dimReductProcessing(dimReductSelect, currTrainingDataStruct);
                
                % then do clustering with circle mapping
                aggregatorSelect.name = 'None';
                obj.classifierModels{ind_motionSet} = ...
                    classifierProcessing(classifierSelect, aggregatorSelect, currTrainingDataStruct, obj.dimReductModels{ind_motionSet});
                obj.classifierModels{ind_motionSet}.adaptiveComponent(currTrainingDataStruct);
            end
            
            obj.motionSet = motionSetEntry;
        end
        
        function [testingData, testingLabel, testingTime] = qdq2stacked(obj, currSubj, dataSettingsSelect, subjectData)
            currLabel = currSubj.label;
            
            permissibleIndStruct = [];
            rawTime = 1:length(currLabel);
            rawLabel = currLabel;
            rawSegTime = [];
            
            objSim.time = 1:length(currLabel);
            objSim.settings = dataSettingsSelect;
            objSim.settings.mode = 'Testing'; % need to create the proper I/Os
            objSim.label = currLabel;
            objSim.subjectData = subjectData;
            rawSegTestingInclude = currSubj.segTestingInclude;
            
            [testingData, testingLabel, ~, ~, testingTime, ~, ~] = ...
                samplingForTestingData(objSim, permissibleIndStruct, rawSegTime, [], rawSegTestingInclude);
        end
        
        function [classifierLabels, groundTruthLabels] = classify(obj,input,extra,blah)
            [dataSettingsSelect, dimReductSelect, classifierSelect] = setSettings(obj);
            
            groundTruthLabels = [];
            classifierLabels = [];
            
            for ind_subjData = 1:length(extra.subjectData)
                % find the proper matrix  
                currSubj = extra.subjectData{ind_subjData};
                currMotionType = currSubj.exerciseType(1:8);
                motionTypeInd = find(strcmpi(currMotionType, obj.motionSet));
                
                % reprocess the data
                [newExtra.testingData, testingLabel, newExtra.testingTime] = ...
                    qdq2stacked(obj, currSubj, dataSettingsSelect, extra.subjectData(ind_subjData));
                groundTruthLabels = [groundTruthLabels; testingLabel];
                
%                 % apply rotation matrix
%                 testingDataPCA = obj.dimReductModels{motionTypeInd}.apply(testingData);
                 
                % apply the classification
                newExtra.subjectData = extra.subjectData(ind_subjData);
                newExtra.useTime = extra.useTime;
                
%                 if isempty(extra.testingTime)
%                     newExtra.testingTime = extra.trainingTime;
%                     newExtra.testingData = extra.trainingData;
%                 else
%                     newExtra.testingTime = extra.testingTime;
%                     newExtra.testingData = extra.testingData;
%                 end
                
%                 localExtra.subjectData = extra.subjectData(ind_subjData);
                testingLabels = obj.classifierModels{motionTypeInd}.classify(newExtra.testingData, newExtra, 0)'; % 0 = no adaptiveness for now
                classifierLabels = [classifierLabels; testingLabels];
            end
        end
        
        function obj_copy = copy(obj)
           obj_copy = SVMClassifier();
           kernels = {'linear','polynomial','radial','sigmoid'};
           obj_copy.init(kernels{obj.kernel+1}); % linear would = 0 if the +1 didn't exist
        end
    end

end