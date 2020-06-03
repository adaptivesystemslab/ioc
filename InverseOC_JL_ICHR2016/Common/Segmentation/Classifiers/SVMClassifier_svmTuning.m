function optimalParam = SVMClassifier_svmTuning(obj, input, output, searchArray)
% iterate though the parameters (first pass)
% all non-RBF kernel parameters are ignored
crossValidationCount = 5;
featureCount = size(input, 2);
degRange = [3];
gammaRange = [1/featureCount searchArray];
coef0Range = [0];
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

function optimalParam = svmTuning(obj, trainingParam, outputDownsample, inputDownsample, crossValidationCount)
% perform parameter search with downsampled data
accuracy = zeros(size(trainingParam));
for ind_trainingParam = 1:length(trainingParam)
    %                  SVMmodel = svmtrain_mex(output,input,['-t ' num2str(obj.kernel) ' -q']);
    
    
    svmConfigStr = [trainingParam(ind_trainingParam).string];
    accuracy(ind_trainingParam) = do_binary_cross_validation(outputDownsample, inputDownsample, svmConfigStr, 5);
    
    %                 svmConfigStr = [trainingParam(ind_trainingParam).string ' -v ' num2str(crossValidationCount)];
    %                 accuracy(ind_trainingParam) = svmtrain_mex(outputDownsample,inputDownsample,svmConfigStr);
    
    fprintf('Tuning SVM (%u/%u)...trying ''%s'': %u', ind_trainingParam, length(trainingParam), svmConfigStr, accuracy(ind_trainingParam));
end

[maxVal, maxInd] = max(accuracy);
optimalParam = trainingParam(maxInd);
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