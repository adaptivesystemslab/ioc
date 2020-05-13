classdef SVMOneClassClassifier < AClassifier
   %This is the Support Vector Machine classifier
   %Uses matlab built in SVM functionality
   
    properties(Access=private)
       SVMmodel = [];
       modelParam = [];
       
       svmprobMode = 'none'; % none platt dist1
       svmprobsigned = 0; % doesn't need to be signed
       svmDecayMode = 0;
       
       svmfwdMode = 0; % 0 no mode, 1 is mode
       svmffsigned = 1; % if 1, then not signed. 
       svmffscaled = 1; % if 1, then not scaled
       
       applyZeroMean = 0; % 1 = remove mean and apply unity covar before training and testing
       applyScaling = 0; % scale the values so the input vector is [0, 1]
       downsampleUpperBound = 10000; % the amount of points used for parameter searching
       meanOffset = [];
       stdDevOffset = [];
       
       %Default Kernel is RBF
       kernel = 2;
       degreeFlag = 0;
       gammaFlag = 0;
       coef0Flag = 0;
       
       ignoreVal = -Inf;
    end
    
    properties(Access=public)
        fsm = [];
    end
    
    methods
        function init(obj,varargin)
            % set up kernel and tuning variables for the SVM training
            
            switch varargin{1}
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
            [input, obj.meanOffset, obj.stdDevOffset] = zeroMeanCovarSignal(input, obj.applyZeroMean, obj.applyZeroMean);
            input = zeroScaleSignal(input, obj.applyScaling);

            % iterate though the parameters (first pass)
%             searchArray = 10.^(-2:2:2);
            searchArray = [];            
            obj.modelParam = determineOptimalParams(obj, input, output, searchArray);
            
             % keep only the p1 points, since one class svm ignores labels
             pointsToKeep = output == 1;
             
             outputToUse = output(pointsToKeep);
             inputToUse = input(pointsToKeep, :);
%              
%              outputToUse = output;
%              inputToUse = input;
%              
%              obj.modelParam.string = '-t 0 -g 0.018182 -c 1 -q -s 2';
%              obj.modelParam.string = '-t 0 -g 0.018182 -c 1 -q -s 0';
%     
            % now train the SVM with all the data
            obj.SVMmodel = svmtrain_mex(outputToUse, inputToUse, obj.modelParam.string);
            
%             obj.plotSupportVectors(inputUnity, output, obj.SVMmodel);
            
            [~, error, ~] =...
                svmpredict_mex(output, input, obj.SVMmodel, '-q');
            error
%             error = 1-error(1)/100;
        end
        
        function trainFsm(obj, trainingData)
          
        end
        
        function optimalParam = determineOptimalParams(obj, input, output, searchArray)
            % iterate though the parameters (first pass)
            crossValidationCount = 5;
            featureCount = size(input, 2);
            degRange = [3 searchArray];
            gammaRange = [1/featureCount searchArray];
            coef0Range = [0 searchArray];
            costRange = [1 searchArray];
            
            % downsample the data for the parameter search if it is too
            % long, to keep the runtime low
            if size(input, 1)/2 > obj.downsampleUpperBound
                % randomly select some points for downsampling
                halfXLength = floor(size(input, 1)/2);
                selectInd = randperm(halfXLength, obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(halfXLength, obj.downsampleUpperBound) + halfXLength;
                inputDownsample = input([selectInd selectInd2], :);
                outputDownsample = output([selectInd selectInd2], :);
                
%                 selectInd = randperm(size(inputUnity, 2), obj.downsampleUpperBound);
%                 inputDownsample = inputUnity(:, selectInd);
            else
                inputDownsample = input;
                outputDownsample = output;
            end
            
            trainingParam = svmTrainString(obj, degRange, gammaRange, coef0Range, costRange);
            
            if ~isempty(searchArray)
                optimalParam = svmTuning(obj, trainingParam, outputDownsample, inputDownsample, crossValidationCount);
                
                % second pass for parameter search
                searchArray = [1];
                degRange = optimalParam.degree*searchArray;
                gammaRange = optimalParam.gamma*searchArray;
                coef0Range =  optimalParam.coef0*searchArray;
                costRange = optimalParam.cost*searchArray;
                
                trainingParam = svmTrainString(obj, degRange, gammaRange, coef0Range, costRange);
                optimalParam = svmTuning(obj, trainingParam, outputDownsample, inputDownsample, crossValidationCount);
            else
                optimalParam = trainingParam;
            end
        end
        
        function plotSupportVectors(obj, input, output, SVMmodel)
            % pull out all the points
            p1 = input(output == 1, :);
            p0 = input(output == 0, :);
            
            % now pull out all the svs
            svIndices = SVMmodel.sv_indices;
            svInput = input(svIndices, :);
            svOutput = output(svIndices);
            svp1 = svInput(svOutput == 1, :);
            svp0 = svInput(svOutput == 0, :);
            
            figure;
            hold on
%             plot(p1(:, 1), p1(:, 2), 'b.');
%             plot(p0(:, 1), p0(:, 2), 'r.');
%             
%             plot(svp1(:, 1), svp1(:, 2), 'bo');
%             plot(svp0(:, 1), svp0(:, 2), 'ro');
            
            skipind = 8;
            p1indtoplot = 1:skipind:size(p1, 1);
            p0indtoplot = 1:skipind:size(p0, 1);
            svp1indtoplot = 1:skipind:size(svp1, 1);
            svp0indtoplot = 1:skipind:size(svp0, 1);
            
%             scatter3(p1(p1indtoplot, 1), p1(p1indtoplot, 2), p1(p1indtoplot, 3), ...
%                 'CData', [0 0 1], 'Marker', '.', 'DisplayName', 'p1 training');
%             
%             scatter3(p0(p0indtoplot, 1), p0(p0indtoplot, 2), p0(p0indtoplot, 3), ...
%                 'CData', [1 0 0], 'Marker', '.', 'DisplayName', 'p0 training');
            
            scatter3(svp1(svp1indtoplot, 1), svp1(svp1indtoplot, 2), svp1(svp1indtoplot, 3), ...
                'CData', [0 0 1], 'Marker', 'o', 'DisplayName', 'p1 sv');
            
            scatter3(svp0(svp0indtoplot, 1), svp0(svp0indtoplot, 2), svp0(svp0indtoplot, 3), ...
                'CData', [1 0 0], 'Marker', 'x', 'DisplayName', 'p0 sv');       
            
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            
            view([-250 10]);
            
            grid on
        end
        
        function trainingParam = svmTrainString(obj, degRange, gammaRange, coef0Range, costRange)
            % based on the obj settings, generate a list of parameters for
            % tuning of the SVM 
            
            % set up the parameter array to iterate through (first value is
            % the 'default' value)
%             degRange = [3 1:2 4:5];
%             gammaRange = [1/featureCount 3/featureCount 5/featureCount];
%             coef0Range = [0 1:5];
%             costRange = [1 2:2:9];
            
            % wide range searches
%             degRange = [3 10.^(-3:3)];
%             gammaRange = [1/featureCount 10.^(-3:3)];
%             coef0Range = [0 10.^(-3:3)];
%             costRange = [1 10.^(-3:3)];
            
            % default
%             degRange = [3];
%             gammaRange = [1/featureCount];
%             coef0Range = [0];
%             costRange = [1];

           if ~obj.degreeFlag 
               degRange = obj.ignoreVal;
           end
           
           if ~obj.gammaFlag 
               gammaRange = obj.ignoreVal;
           end
           
           if ~obj.coef0Flag
               coef0Range = obj.ignoreVal;
           end

            trainingParamCounter = 0;
            for ind_degree = 1:length(degRange)
                for ind_gamma = 1:length(gammaRange)
                    for ind_coef0 = 1:length(coef0Range)
                        for ind_cost = 1:length(costRange)
                            trainingParamInstance.kernel = obj.kernel;
                            trainingParamInstance.degree = degRange(ind_degree);
                            trainingParamInstance.gamma = gammaRange(ind_gamma);
                            trainingParamInstance.coef0 = coef0Range(ind_coef0);
                            
                            trainingParamInstance.cost = costRange(ind_cost);
                            
                            trainingParamInstance.string = convertToStringParam(obj, trainingParamInstance);
            
                            trainingParamCounter = trainingParamCounter + 1;
                            trainingParam(trainingParamCounter) = trainingParamInstance;
                        end
                    end
                end
            end
        end
        
        function [out, prob] = classify(obj,input,labels,incMode,extra)
            input = preprocessData(obj, input);
            
            % shift the data to 0 mean, 1 covar if it is set as a flag
            input = zeroMeanCovarSignal(input, obj.applyZeroMean, obj.applyZeroMean, obj.meanOffset, obj.stdDevOffset);
            input = zeroScaleSignal(input, obj.applyScaling);
            
            [out, ~, decision_values] = svmpredict_mex(zeros(size(input,1),1), input, obj.SVMmodel, '-q');
            
            out(out == -1) = 0;
            out(out ==  1) = 1;
%             
%             out(out ==  1) = 0;
%             out(out == -1) = 1;
                        
            
            prob = ones(size(out))*0.5;
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
        
        function instanceString = convertToStringParam(obj, trainingParamInstance)
            % define the kernel
            instanceString = ['-t ' num2str(trainingParamInstance.kernel)];
            
            if ~isequal(trainingParamInstance.degree, obj.ignoreVal)
                instanceString = [instanceString ' -d ' num2str(trainingParamInstance.degree)];
            end
            
            if ~isequal(trainingParamInstance.gamma, obj.ignoreVal)
                instanceString = [instanceString ' -g ' num2str(trainingParamInstance.gamma)];
            end
            
            if ~isequal(trainingParamInstance.coef0, obj.ignoreVal)
                instanceString = [instanceString ' -r ' num2str(trainingParamInstance.coef0)];
            end
            
            if ~isequal(trainingParamInstance.cost, obj.ignoreVal)
                instanceString = [instanceString ' -c ' num2str(trainingParamInstance.cost)];
            end
            
            % then add quiet mode
            instanceString = [instanceString ' -q'];
            
            % then add single class mode
%             instanceString = [instanceString ' -s 0']; 
            instanceString = [instanceString ' -s 2']; 
            
            % ['-t ' num2str(obj.kernel) ' -q']
        end
        
        function optimalParam = svmTuning(obj, trainingParam, outputDownsample, inputDownsample, crossValidationCount)
            % perform parameter search with downsampled data
            accuracy = zeros(size(trainingParam));
            for ind_trainingParam = 1:length(trainingParam)
%                  SVMmodel = svmtrain_mex(output,input,['-t ' num2str(obj.kernel) ' -q']);
                svmConfigStr = [trainingParam(ind_trainingParam).string ' -v ' num2str(crossValidationCount)];
                fprintf('Tuning SVM (%u/%u)...trying ''%s'': ', ind_trainingParam, length(trainingParam), svmConfigStr);
                accuracy(ind_trainingParam) = svmtrain_mex(outputDownsample,inputDownsample,svmConfigStr);
            end
            
            [maxVal, maxInd] = max(accuracy);
            optimalParam = trainingParam(maxInd);
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