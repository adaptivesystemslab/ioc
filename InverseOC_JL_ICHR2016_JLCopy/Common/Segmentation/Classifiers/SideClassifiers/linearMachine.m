classdef linearMachine < AClassifier
% FSM
    properties(SetAccess=protected)
        k = 2;
        dataFormat = []; % which data to use?
        
        trainingData = []; % training data
        trainingLabel = []; % the labels 
        trainingLabelGndTruth = []; % the ground truth
        
        trainingDataOrig = [];
        trainingLabelOrig = [];
        trainingLabelGndTruthOrig = [];
        
        outputUnique = [];
        clusterCentreArray = [];
        clusterDistanceArray = [];
        clusterLabelArray = [];
        clusterCovarArray = [];
        
        dimReduct = [];
        
        training = [];
        group = [];
        
        fsm = [];
        
        lastTrainedTimeStamp = [];
        
        %     auxState.rawRest = []; % 1 resting
        %     auxState.rawFlex = []; % 2 moving towards peak
        %     auxState.rawPeak = []; % 3 around peak
        %     auxState.rawExt = [];  % 4 coming back from peak
        labels = [1 0 1 0];
    end
    
    methods
        function obj = linearMachine(incParam)
            obj.dataFormat = incParam.dataFormat;
            obj.k = incParam.clusterCountPerClass;
        end
            
        function init(obj,varargin)
        % Initialization function for FSM

        end
        
        function updateDimReductClass(obj, dimReduct)
            % save the dim reduction into the class
            obj.dimReduct = dimReduct;
        end
        
        function [clusterCentre, clusterDistance, clusterCovar] = calcKmeans(obj, data, k)
            % calculate kmeans
            [IDX, clusterCentre, SUMD, D] = kmeans(data, k);
            
            dist = sqrt(D); % kmeans distance is reported to be the sq euclidean
            
%             for i = 1:k
%                 distVec(:, i) = normVector(repmat(clusterCentre(i, :), size(data, 1), 1) - data);
%             end
            
            % determine distance
            for i = 1:k
                distK = dist(IDX == i, i);
                clusterDistance(i) = max(distK);
                
                entryK = data(IDX == i, :);
                clusterCovar{i} = cov(entryK);
            end
        end
        
        function [labelArrayOut, labelArrayProb] = labelOnDistance(obj, fullDataInput)
            % calculate the distance to each point
            for i = 1:size(obj.clusterCentreArray, 1)
                currCentre = obj.clusterCentreArray(i, :); 
                currCovar = obj.clusterCovarArray{i};
                
%                 for j = 1:size(fullDataInput, 1)
%                     dataInput = fullDataInput(j, :);
%                     dataInputCentred = dataInput - repmat(currCentre, size(dataInput, 1), 1);
%                     dataInputCentred = dataInputCentred';
%                     dataInputDist(j) = dataInputCentred' * inv(currCovar) * dataInputCentred;
%                 end
            
                dataInput = fullDataInput(:, :);
                dataInputCentred = dataInput - repmat(currCentre, size(dataInput, 1), 1);
                dataInputCentred = dataInputCentred';
                dataInputDist2 = dataInputCentred' * inv(currCovar) * dataInputCentred;
                dataInputDist = diag(dataInputDist2);

                labelArray(:, i) = dataInputDist;
            end
            
            for i = 1:size(labelArray)
                [labelDist, labelInd] = min(labelArray(i, :));
                labelArrayOut(i) = obj.clusterLabelArray(labelInd);
                labelArrayProb(i) = 0.5;
            end
            
%             % go for the smallest distance
%             for i = 1:size(labelArray)
%                 currLabelArray = labelArray(i, :);
%                 [C1, I1] = min(currLabelArray);
%                 currLabelArray(I1) = currLabelArray(I1) + 10000;
%                 [C2, I2] = min(currLabelArray);
%                 
%                 if C1 > 0.85
%                     % if the min is more than 85% of the way away, ie we're uncertain
%                     labelToApply = -1;
%                     probToApply = 1;
%                 elseif (obj.clusterLabelArray(I1) ~= obj.clusterLabelArray(I2) ...
%                         && GndTruth(C2 - C1) < 0.10) 
%                     % two different labels, and are close
%                     labelToApply = -1;
%                     probToApply = 1;
%                 else
%                     labelToApply = obj.clusterLabelArray(I1);
%                     probToApply = 1;
%                 end
%                 
%                 labelArrayOut(i) = labelToApply;
%                 labelArrayProb(i) = C1;
%             end
            
%             % determine the label. if it's in an overlapping region, assign
%             % it an unknown label
%             for i = 1:length(obj.outputUnique)
%                 labelArrayCombined(:, i) = sum(labelArray(:, obj.clusterLabelArray == obj.outputUnique(i)), 2);
%             end
%             
%             for i = 1:size(labelArrayCombined)
%                 currLabelArrayTest = find(labelArrayCombined);
%                 if length(currLabelArrayTest == 1)
%                     labelArrayOut(i) = labelArrayCombined(currLabelArrayTest);
%                 else
%                     labelArrayOut(i) = -1;
%                 end
%             end
        end
        
        function error = retrain(obj,dataStructInput,percentageToKeep,distanceToKeep)
            % combine the original input with the new output
            
            numberOfData = floor(length(obj.trainingLabel) * percentageToKeep);
            numberOfDataToKeep = randperm(length(obj.trainingLabel), numberOfData);
            
            if ~isempty(dataStructInput)
                newInput = [obj.trainingData(numberOfDataToKeep, :); dataStructInput.input];
                newOutput = [obj.trainingLabel(numberOfDataToKeep, :); dataStructInput.output];
                newOutputGndTruth = [obj.trainingLabelGndTruth(numberOfDataToKeep, :); dataStructInput.outputGndTruth];
            else
                % if passed in a blank, ignore the state
                newInput = obj.trainingData(numberOfDataToKeep, :);
                newOutput = obj.trainingLabel(numberOfDataToKeep, :);
                newOutputGndTruth = obj.trainingLabelGndTruth(numberOfDataToKeep, :);
            end
            
            dataStruct.input = newInput; % the original data
            dataStruct.output = newOutput;
            dataStruct.outputGndTruth = newOutputGndTruth;
            
            train(obj, dataStruct, distanceToKeep);
        end
        
        function resetTrainingData(obj)
%             train(obj, obj.trainingDataOrig, obj.trainingLabelOrig);
            
            obj.trainingData = obj.trainingDataOrig;
            obj.trainingLabel = obj.trainingLabelOrig;
            obj.trainingLabelGndTruth = obj.trainingLabelGndTruthOrig;
        end
        
        function inputDR = applyDataFormat(obj, input, output)
           switch obj.dataFormat
                case 'PCA'
                    % now perform the dim reduction, since the DR may change from
                    % instance to instance
                    inputDR = obj.dimReduct.apply(input);
                    inputDR = inputDR(:, 1:2);
                    
                case 'PP'                    
                    dofs = 1:size(input, 2)/2;
                    
                    xdof = normVector(input(:, dofs));
                    ydof = normVector(input(:, dofs+max(dofs)));
                    
                    inputDR = [xdof ydof];
           end 
            
%            x_shifted = inputDR(:, 1);
%            y_shifted = inputDR(:, 2);
                    x_shifted = inputDR(:, 1) - mean(inputDR(:, 1));
                    y_shifted = inputDR(:, 2) - mean(inputDR(:, 2));
                    
                    x_shifted = abs(x_shifted);
                    y_shifted = abs(y_shifted);
                    
                    [polarTheta,polarRho] = cart2pol(x_shifted,y_shifted);
                    
                    if 0
                        t = 1:length(output); 
                        h1 = figure('position', [       1384         147        1240         913]);
                        
                        h1(1) = subplot(2, 1, 1);
                        plot(x_shifted(output == 1), y_shifted(output == 1), 'r.'); hold on
                        plot(x_shifted(output == 0), y_shifted(output == 0), 'b.');
                        h1(2) = subplot(2, 1, 2);
                        polar(polarTheta(output == 1), polarRho(output == 1), 'r.'); hold on
                        polar(polarTheta(output == 0), polarRho(output == 0), 'b.')
                        
%                         h1(1) = subplot(2, 1, 1);
% %                       plot(t(output == 1), input(output == 1, dofs), 'r.'); hold on
% %                         plot(t(output == 0), input(output == 0, dofs), 'b.');
%                         
%  h1(2) = subplot(2, 1, 2);
% %                         plot(t(output == 1), polarTheta(output == 1), 'r.'); hold on
%                         plot(t(output == 0), polarTheta(output == 0), 'b.');
%                                       linkaxes(h1, 'x');
          
                        
                        close(h1)
                    end
                    
                    inputDR = polarTheta;
        end
        
        function error = train(obj,dataStruct,distanceToKeep)
            if ~exist('percentToKeep', 'var')
                % distance from core to keep
                distanceToKeep = 1;
            end
            
            input = dataStruct.input;
            output = dataStruct.output;
            outputGndTruth = dataStruct.outputGndTruth;
            
            inputDR = applyDataFormat(obj, input, output);
            
            if isempty(obj.lastTrainedTimeStamp)
                % first time training. copy out the original data in case
                % we need to reset the classifier. make sure to save it in
                % the original space and not DR
                obj.trainingDataOrig = input;
                obj.trainingLabelOrig = output;
                obj.trainingLabelGndTruthOrig = outputGndTruth;
            end

            
            % pull out the data for the two classes
            obj.outputUnique = unique(output);
            
            for i = 1:length(obj.outputUnique)
                outputSelection = output == obj.outputUnique(i);
                numberOfOutput = length(find(outputSelection));
                dataUse{i} = inputDR(outputSelection, :);
                
                if ~isempty(obj.k) && numberOfOutput > 1
                    % calculate mean and distance
                    [clusterCentre{i}, clusterDistance{i}, clusterCovar{i}] = calcKmeans(obj, dataUse{i}, obj.k);
%                 elseif numberOfOutput == 1
%                     clusterCentre{i} = dataUse{i};
%                     clusterDistance{i} = 1;
                else
                    % if there's only one input, we can't normalize for
                    % distance, so we'll 
                    clusterCentre{i} = [];
                    clusterCovar{i} = [];
                    clusterDistance{i} = [];
                end
                
                obj.clusterCentreArray = [obj.clusterCentreArray; clusterCentre{i}];
                obj.clusterDistanceArray = [obj.clusterDistanceArray; clusterDistance{i}'];
                obj.clusterCovarArray = [obj.clusterCovarArray; clusterCovar{i}'];
                obj.clusterLabelArray = [obj.clusterLabelArray; (ones(size(clusterDistance{i}))*obj.outputUnique(i))'];
            end

            % test
            [labelArrayOut, labelArrayProb] = classify(obj, input);
            err = sum(labelArrayOut' == output)/length(output);
            
            if 0
                labelArrayOut = labelArrayOut';
                h = figure('position', [  1375         303         560*2         420]);
                color = {'r', 'g', 'b'};
                t = [1:length(labelArrayOut)]';                
                
                for i = 1:length(obj.outputUnique)
                    subplot(1, 2, 1);
                    hold on
                    lab = obj.outputUnique(i);
                    labInd1 = find(labelArrayOut == lab);
                    plot(inputDR(labInd1, 1), inputDR(labInd1, 2), [color{i} '.'], 'DisplayName', 'ClusterAssignment');
%                     plot(t(labInd1), inputDR(labInd1), [color{i} '.'], 'DisplayName', 'ClusterAssignment');
                    title('Cluster');
                    
                    subplot(1, 2, 2);
                    hold on
                    labInd2 = find(outputGndTruth == lab);
                    plot(inputDR(labInd2, 1), inputDR(labInd2, 2), [color{i} '.'], 'DisplayName', 'GroundTruth');
%                     plot(t(labInd2), inputDR(labInd2), [color{i} '.'], 'DisplayName', 'GroundTruth');
                    title('Ground');
                end
                
%                 legend('show')
                
                for i = 1:length(obj.outputUnique)
                    subplot(1, 2, 1);  hold on
                    plot(clusterCentre{i}(:, 1), clusterCentre{i}(:, 2), 'kx', 'markersize', 20);
                    error_ellipse(clusterCovar{i}{1}(1:2, 1:2), clusterCentre{i}(1, 1:2));
                    error_ellipse(clusterCovar{i}{2}(1:2, 1:2), clusterCentre{i}(2, 1:2));

%                     plot(t, ones(size(t))*clusterCentre{i}(1), 'kx');
%                     plot(t, ones(size(t))*clusterCentre{i}(2), 'kx');
                    
                    subplot(1, 2, 2);  hold on
                    plot(clusterCentre{i}(:, 1), clusterCentre{i}(:, 2), 'kx', 'markersize', 20);
                    error_ellipse(clusterCovar{i}{1}(1:2, 1:2), clusterCentre{i}(1, 1:2));
                    error_ellipse(clusterCovar{i}{2}(1:2, 1:2), clusterCentre{i}(2, 1:2));

%                     plot(t, ones(size(t))*clusterCentre{i}(1), 'kx');
%                     plot(t, ones(size(t))*clusterCentre{i}(2), 'kx');
                end
                
                subplot(1, 2, 1); xlim([-6 6]); ylim([-6 6]);
                subplot(1, 2, 2); xlim([-6 6]); ylim([-6 6]); 
                
                saveas(h, ['C:\Documents\MATLABResults\imgDump\' datestr(now, 'yyyy-mm-dd-HH-MM-SS') '.jpg']);
%                 close(h);
            end
            
            % prune the dataset in prep for future runs
            % if the labelArrayProb says it's within n% of a cluster core,
            % let's keep it
            probToKeep = labelArrayProb <= distanceToKeep;
            obj.trainingData = input(probToKeep, :); % save the data in the ORIG space, not DR space
            obj.trainingLabel = output(probToKeep);
            obj.trainingLabelGndTruth = outputGndTruth(probToKeep);
            
            obj.lastTrainedTimeStamp = datestr(now, 'yyyy-mm-dd HH-MM-SS');
        end
        
        function [out, prob] = classify(obj,input)
            jumps = 1000;
            for i = 1:jumps:size(input, 1) % memory overflow issues when the array is too big
                indStart = i;
                indEnd = i+jumps;
                if indEnd > size(input, 1)
                    indEnd = size(input, 1);
                end
                
                inputDR = applyDataFormat(obj, input(indStart:indEnd, :));
                [out(indStart:indEnd), prob(indStart:indEnd)] = labelOnDistance(obj, inputDR);
            end
        end
        
        function obj_copy = copy(obj)
            obj_copy = finiteStateMachine();
        end
    end    
end