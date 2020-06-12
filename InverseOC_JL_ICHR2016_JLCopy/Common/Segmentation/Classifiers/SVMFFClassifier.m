classdef SVMFFClassifier < AClassifier
   %This is the Support Vector Machine classifier
   %Uses matlab built in SVM functionality
   
    properties(Access=private)
       SVMmodel = [];
       modelParam = [];
       
       svmprobMode = 'none'; % none platt dist1
       svmfwdMode = 0; % 0 no mode, 1 is mode
       
       downsampleUpperBound = 1500; % the amount of points used for parameter searching
       applyZeroMean = 0; % 1 = remove mean and apply unity covar before training and testing
       meanOffset = [];
       stdDevOffset = [];
       
       %Default Kernel is RBF
       kernel = 2;
       degreeFlag = 0;
       gammaFlag = 0;
       coef0Flag = 0;
       
       ignoreVal = -Inf;
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
                [inputUnity, obj.meanOffset, obj.stdDevOffset] = zeroMeanCovarSignal(input, 1, 1);
            else
                inputUnity = input; 
            end
            
            % downsample the data for the parameter search if it is too
            % long, to keep the runtime low
            if size(inputUnity, 1)/2 > obj.downsampleUpperBound
                % randomly select some points for downsampling
                halfXLength = floor(size(inputUnity, 1)/2);
                selectInd = randperm(halfXLength, obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(halfXLength, obj.downsampleUpperBound) + halfXLength;
                inputDownsample = inputUnity([selectInd selectInd2], :);
                outputDownsample = output([selectInd selectInd2], :);
                
%                 selectInd = randperm(size(inputUnity, 2), obj.downsampleUpperBound);
%                 inputDownsample = inputUnity(:, selectInd);
            else
                inputDownsample = inputUnity;
                outputDownsample = output;
            end
            
            if obj.svmfwdMode
                % augment the features with the forward feeding mode
                outputLog = ones(size(outputDownsample)); % all values that cannot be assessed are labelled 1
                counter = 1;
                
                for ind_input = 1:length(outputLog)
                    if ind_input == 1
                        counter = 1;
                    elseif ind_input == 2
                        counter = 2;
                    else
                        counter = labelChecker(obj, counter, output(ind_input-1), output(ind_input-2));
                    end
                    
                    if ind_input > 1 && output(ind_input-1) < 1
                        outputLog(ind_input) = -counter/100;
                    else
                        outputLog(ind_input) = counter/100;
                    end
                end
                
                inputDownsample = [inputDownsample outputLog];
            end
            
            % iterate though the parameters (first pass)
            crossValidationCount = 5;
            featureCount = size(input, 2);
%             searchArray = 10.^(-2:2:2);
            searchArray = [];
            degRange = [3 searchArray];
            gammaRange = [1/featureCount searchArray];
            coef0Range = [0 searchArray];
            costRange = [1 searchArray];
            
            trainingParam = svmTrainString(obj, degRange, gammaRange, coef0Range, costRange);
            optimalParam = svmTuning(obj, trainingParam, outputDownsample, inputDownsample, crossValidationCount);
            
            % second pass for parameter search
            searchArray = [1];
            degRange = optimalParam.degree*searchArray;
            gammaRange = optimalParam.gamma*searchArray;
            coef0Range =  optimalParam.coef0*searchArray;
            costRange = optimalParam.cost*searchArray;
            
            trainingParam = svmTrainString(obj, degRange, gammaRange, coef0Range, costRange);
            optimalParam = svmTuning(obj, trainingParam, outputDownsample, inputDownsample, crossValidationCount);
            
            obj.modelParam = optimalParam;
            
            switch obj.svmprobMode
                case 'none'
                    
                case 'platt'
                    obj.modelParam.string = [obj.modelParam.string ' -b 1'];
                case 'dist1'
                    
            end
   
            % now train the SVM with the zero mean data
            obj.modelParam
            obj.SVMmodel = svmtrain_mex(output,inputUnity,obj.modelParam.string);
            
%             obj.plotSupportVectors(inputUnity, output, obj.SVMmodel);
            
            [~, error, ~] =...
                svmpredict_mex(output, input, obj.SVMmodel, '-q');
            error = 1-error(1)/100;
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
        
        function [out, prob] = classify(obj,input,labels,extra)
            input = preprocessData(obj, input);
            
            if obj.applyZeroMean
                % shift the data to 0 mean, 1 covar
                inputUnity = zeroMeanCovarSignal(input, 1, 1, obj.meanOffset, obj.stdDevOffset);
            else
                inputUnity = input; 
            end
            
            if obj.svmfwdMode
                out = -1*ones(size(input, 1), 1);
                outputLog = out;
                decision_values = -1*ones(size(input, 1), 1);
                
                for ind_input = 1:size(input, 1)
                    if ind_input == 1
                        counter = 1;
                    elseif ind_input == 2
                        counter = 2;
                    else
                        counter = labelChecker(obj, counter, out(ind_input-1), out(ind_input-2));
                    end
                    
                    if ind_input > 1 && out(ind_input-1) < 1
                        outputLog(ind_input) = -counter/100;
                    else
                        outputLog(ind_input) = counter/100;
                    end
                    
                    currentInput = [input(ind_input, :) outputLog(ind_input)];
                    
                    [out(ind_input), ~, decision_values(ind_input)] = ...
                        svmpredict_mex(zeros(size(currentInput,1),1), currentInput, obj.SVMmodel,'-q');
                end
            else
                [out, ~, decision_values] = svmpredict_mex(zeros(size(inputUnity,1),1), inputUnity, obj.SVMmodel,'-q');
            end
            
            switch obj.svmprobMode
                case 'none'
                    prob = ones(size(decision_values));
                    
                case 'platt'
                    dist = decision_values;
                    prob = sigmoid_predict(obj, decision_values, obj.SVMmodel.ProbA, obj.SVMmodel.ProbB);
                case 'dist1'
                    % calculating distance of point to boundary
                    % http://stats.stackexchange.com/questions/126379/libsvm-on-matlab-with-rbf-kernel-compute-distance-from-hyperplane
                    gamma = obj.modelParam.gamma;
                    rho = obj.SVMmodel.rho;
                    SV = obj.SVMmodel.SVs;
                    sv_coef = obj.SVMmodel.sv_coef;
                    
                    for i = 1:size(inputUnity, 1)
                        currObs = inputUnity(i, :);
                        for k = 1:size(SV, 1)
                            expVal = SV(k) - currObs;
                            normCoeff = norm(expVal);
                            kernelize(k) = sv_coef(k) * exp(-gamma * normCoeff^2);
                        end
                        dist(i) = sum(kernelize) - rho;
                    end
                    
%                     dist2 = decision_values ./ dist';
                    
                    dist(out == 0) = -dist(out == 0);
                    
                    % running it through a sigmoid
                    A = 5;
                    B = 0;
                    
%                     dist = -10:0.1:10;
                    prob = 1./(1 + exp(-A*(dist-B))) - 0.5;
%                     plot(dist, sig);
            end
            
            
% % %             if ~strcmpi(obj.svmprobMode, 'none')
% % %                 % if there is some probability, use the probability to
% % %                 % modify the outputs
% % %                 probWindowLen = 4;
% % %                 weightedProbScale = [0.1 0.2 0.3 0.4];
% % %                 probGapThreshold = 0.5;
% % %                 
% % %                 prob_new = prob;
% % %                 out_new = out;
% % %                 
% % %                 for ind_labels = probWindowLen+1:length(prob)
% % %                     probWindow = prob(ind_labels-probWindowLen:ind_labels-1);
% % %                     labelWindow = out(ind_labels-probWindowLen:ind_labels-1);
% % %                     probCurr = prob(ind_labels);
% % %                     labelCurr = out(ind_labels);
% % %                     
% % %                     if abs(probCurr) < 1 ... % if ambigious to what the label is, and current
% % %                             && labelWindow(end) ~= labelCurr ... % label doesn't match up with previous label
% % %                             && abs(probCurr - probWindow(end)) < probGapThreshold % and the prob jump is small, emphasizing ambigious labels
% % %                         
% % %                         % if the weighted prob of the previous labels is
% % %                         % higher than the current one, then maybe it's
% % %                         % worth flipping the label
% % %                         weightedProb = sum(abs(weightedProbScale * probWindow));
% % %                         if weightedProb > abs(probCurr)
% % %                             % flip the label
% % %                             
% % %                              if out(ind_labels) == 0
% % %                                  out_new(ind_labels) = 1;
% % %                              else
% % %                                  out_new(ind_labels) = 0;
% % %                             end
% % %                             
% % %                             prob_new(ind_labels) = weightedProb * sign(prob(ind_labels));
% % %                         end
% % %                     end
% % %                 end
% % %                 % 
% % %             end
           
            if 0
                if exist('extra', 'var')
                    testingTime = extra.testingTime - extra.testingTime(1);
                    startTime = extra.testingSegTime.startTimeVal - extra.testingTime(1);
                    endTime = extra.testingSegTime.endTimeVal - extra.testingTime(1);
                    titleStr = extra.identString;
                else
                    testingTime = 1:length(prob);
                    startTime = [];
                    endTime = [];
                    titleStr = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
                end
                
                % plot the stuff
            h =    figure;
                axes(1) = subplot(311);
                plot(testingTime, inputUnity(:, 1:5));
                plotBoxes(h, startTime, endTime, 'r', 0, 3, -3);
                        
%                 hold on
%                 plot(out, 'bx');
%                 if ~strcmpi(obj.svmprobMode, 'none')
% %                     plot(prob_new, 'r');
%                     plot(out_new, 'rx');
%                 end
                
                title(titleStr);
                
                axes(2) = subplot(312);
                plot(testingTime, dist, 'b');
                hold on
                plot(testingTime, out, 'bx');
                if exist('labels', 'var')
                    plot(testingTime, labels, 'ro');
                end
                title('dist');          
                
%                 axes(2) = subplot(312);
%                 plot(testingTime, prob, 'b');
%                 hold on
%                 plot(testingTime, out, 'bx');
%                 if exist('labels', 'var')
%                       plot(testingTime, labels, 'ro');
%                 end
%                 title('prob (old) - classifier (b) and ground truth (r)');
                
                axes(3) = subplot(313);
                plot(testingTime, prob, 'b');
                hold on
                plot(testingTime, out, 'bx');
                title('prob (new) - classifier (b) and ground truth (r)');
%                 title('prob (new)');
                if exist('labels', 'var')
                    plot(testingTime, labels, 'ro');
                end
                ylim([-0.1 1.1]);
                linkaxes(axes, 'x');
                
                saveas(h, ['D:\MATLABResults\Seg_Healthy1\figs\' titleStr '.fig']);
                close(h);
            end
        end
        
        function emission = sigmoid_predict(obj, decision_value, A, B)
            fApB = decision_value*A+B;
            
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