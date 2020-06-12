classdef SVMClassifier < AClassifier
   %This is the Support Vector Machine classifier
   %Uses matlab built in SVM functionality
   
    properties(Access=public)
       SVMmodel = [];
       modelParam = [];
       
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
       
       searchArray = 10.^(-3:1:3);
%        searchArray = [];
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
            
            switch classifierSelect.tuning % tuning
                case '1'
                    obj.searchArray = 10.^(-3:1:3);
                case '0'
                    obj.searchArray = [];
            end
        end
        
        function initializeAdaptiveElements(obj)
            
        end
        
        function input = preprocessData(obj, input)
%             % exp: convert to cylindrical
%             [THETA,RHO,Z] = cart2pol(input(:, 1), input(:, 2), input(:, 3));
%             input = [THETA, RHO, Z];

            % convert to polar coordinates
            
        end
        
        function error = train(obj,input,output)
            % this function performs a simple parameter search, then trains
            % the SVM based on the optimal parameter
            input = preprocessData(obj, input);
            
            % for SVM, the data shold be normalized before inserting into
            % the system: zero mean and unity covariance. store the offset
            % covar to apply to the testing data
            if obj.applyZeroMean
                [input, obj.meanOffset, obj.stdDevOffset] = zeroMeanCovarSignal(input, 1, 1);
            end
            
            % iterate though the parameters (first pass)           
            obj.modelParam = SVMClassifier_svmTuning(obj, input, output, obj.searchArray);
            
            % now train the SVM with all the data
            obj.SVMmodel = svmtrain_mex(output, input, obj.modelParam.string);
            
            [trainingAlgLabel, error, ~] = svmpredict_mex(output, input, obj.SVMmodel, '-q');
            error = 1-error(1)/100;
            acc = sum(output == trainingAlgLabel)/length(trainingAlgLabel);
            %             obj.plotSupportVectors(inputUnity, output, obj.SVMmodel);
        end
        
        function [out, prob] = classify(obj,input,labels,incMode,extra)
            input = preprocessData(obj, input);
            
            if obj.applyZeroMean
                % shift the data to 0 mean, 1 covar
                input = zeroMeanCovarSignal(input, 1, 1, obj.meanOffset, obj.stdDevOffset); 
            end
            
            [out, ~, decision_values] = svmpredict_mex(zeros(size(input,1),1), input, obj.SVMmodel, '-q');
            prob = ones(size(decision_values)) * 0.5;         
        end
        
        function obj_copy = copy(obj)
           obj_copy = SVMClassifier();
           kernels = {'linear','polynomial','radial','sigmoid'};
           obj_copy.init(kernels{obj.kernel+1}); % linear would = 0 if the +1 didn't exist
        end
    end

end