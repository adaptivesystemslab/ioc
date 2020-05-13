classdef clusterMachine < AClassifier
% FSM
    properties(SetAccess=protected)
        k = 2;
        
        trainingData = [];
        trainingLabel = [];        
        
        trainingDataOrig = [];
        trainingLabelOrig = [];
        
        outputUnique = [];
        clusterCentreArray = [];
        clusterDistanceArray = [];
        clusterLabelArray = [];
        
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
        function obj = clusterMachine(incParam)
            obj.k = incParam.clusterCountPerClass;
        end
            
        function init(obj,varargin)
        % Initialization function for FSM

        end
        
        function updateDimReductClass(obj, dimReduct)
            % save the dim reduction into the class
            obj.dimReduct = dimReduct;
        end
        
        function [clusterCentre, clusterDistance] = calcKmeans(obj, data, k)
            % calculate kmeans
            [IDX, clusterCentre, SUMD, D] = kmeans(data, k);
            
            dist = sqrt(D); % kmeans distance is reported to be the sq euclidean
            
%             for i = 1:k
%                 distVec(:, i) = normVector(repmat(clusterCentre(i, :), size(data, 1), 1) - data);
%             end
            
            % determine distance
            for i = 1:k
                clusterDistance(i) = max(dist(IDX == i, i));
            end
        end
        
        function [labelArrayOut, labelArrayProb] = labelOnDistance(obj, dataInput)
            % calculate the distance to each point
            for i = 1:size(obj.clusterCentreArray, 1)
                currCluster = obj.clusterCentreArray(i, :);
                distVec = repmat(currCluster, size(dataInput, 1), 1) - dataInput;
                dist = normVector(distVec);
                
                % normalize against the max distance. the smaller the
                % number, the more closer it is
                currLabelArray = 1 - (obj.clusterDistanceArray(i) - dist)/obj.clusterDistanceArray(i); 
                labelArray(:, i) = currLabelArray;
            end
            
            % go for the smallest distance
            for i = 1:size(labelArray)
                currLabelArray = labelArray(i, :);
                [C1, I1] = min(currLabelArray);
                currLabelArray(I1) = currLabelArray(I1) + 100;
                [C2, I2] = min(currLabelArray);
                
                if C1 > 0.85
                    % if the min is more than 85% of the way away, ie we're uncertain
                    labelToApply = -1;
                    probToApply = 1;
                elseif (obj.clusterLabelArray(I1) ~= obj.clusterLabelArray(I2) ...
                        && abs(C2 - C1) < 0.10) 
                    % two different labels, and are close
                    labelToApply = -1;
                    probToApply = 1;
                else
                    labelToApply = obj.clusterLabelArray(I1);
                    probToApply = 1;
                end
                
                labelArrayOut(i) = labelToApply;
                labelArrayProb(i) = C1;
            end
            
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
        
        function error = retrain(obj, input,output, percentageToKeep, distanceToKeep)
            % combine the original input with the new output
            
            numberOfData = floor(length(obj.trainingLabel) * percentageToKeep);
            numberOfDataToKeep = randperm(length(obj.trainingLabel), numberOfData);
            
            newInput = [input; obj.trainingData(numberOfDataToKeep, :)];
            newOutput = [output; obj.trainingLabel(numberOfDataToKeep, :)];
            
            train(obj, newInput, newOutput, distanceToKeep);
        end
        
        function resetTrainingData(obj)
%             train(obj, obj.trainingDataOrig, obj.trainingLabelOrig);
            
            obj.trainingData = obj.trainingDataOrig;
            obj.trainingLabel = obj.trainingLabelOrig;
        end
        
        function error = train(obj,input,output,distanceToKeep)
            
            if ~exist('percentToKeep', 'var')
                % distance from core to keep
                distanceToKeep = 0.9;
            end
            
            if isempty(obj.lastTrainedTimeStamp)
                % first time training. copy out the original data in case
                % we need to reset the classifier. make sure to save it in
                % the original space and not DR
                obj.trainingDataOrig = input;
                obj.trainingLabelOrig = output;
            end
            
            % now perform the dim reduction, since the DR may change from
            % instance to instance
            inputDR = obj.dimReduct.apply(input);
            
            % pull out the data for the two classes
            obj.outputUnique = unique(output);
            
            for i = 1:length(obj.outputUnique)
                outputSelection = output == obj.outputUnique(i);
                numberOfOutput = length(find(outputSelection));
                dataUse{i} = inputDR(outputSelection, :);
                
                if ~isempty(obj.k) && numberOfOutput > 1
                    % calculate mean and distance
                    [clusterCentre{i}, clusterDistance{i}] = calcKmeans(obj, dataUse{i}, obj.k);
%                 elseif numberOfOutput == 1
%                     clusterCentre{i} = dataUse{i};
%                     clusterDistance{i} = 1;
                else
                    % if there's only one input, we can't normalize for
                    % distance, so we'll 
                    clusterCentre{i} = [];
                    clusterDistance{i} = [];
                end
                
                obj.clusterCentreArray = [obj.clusterCentreArray; clusterCentre{i}];
                obj.clusterDistanceArray = [obj.clusterDistanceArray; clusterDistance{i}'];
                obj.clusterLabelArray = [obj.clusterLabelArray; (ones(size(clusterDistance{i}))*obj.outputUnique(i))'];
            end
            
         
%             % cluster
%             figure;
%             color = {'r', 'g', 'b'};
%             for i = 1:length(obj.outputUnique)
%                 hold on
%                 plot(dataUse{i}(:, 1), dataUse{i}(:, 2), [color{i} '.']);
%                 plot(clusterCentre{i}(:, 1), clusterCentre{i}(:, 2), 'kx', 'markersize', 15);
%             end
%             xlim([-10 10]); ylim([-8 8]);
            
            % test
            [labelArrayOut, labelArrayProb] = classify(obj, input);
            err = sum(labelArrayOut' == output)/length(output);
            
            % prune the dataset in prep for future runs
            % if the labelArrayProb says it's within n% of a cluster core,
            % let's keep it
            probToKeep = labelArrayProb < distanceToKeep;
            obj.trainingData = input(probToKeep, :); % save the data in the ORIG space, not DR space
            obj.trainingLabel = output(probToKeep);
            
            obj.lastTrainedTimeStamp = datestr(now, 'yyyy-mm-dd HH-MM-SS');
        end
        
        function [out, prob] = classify(obj,input)
            inputDR = obj.dimReduct.apply(input);
            [out, prob] = labelOnDistance(obj, inputDR);
        end
        
        function obj_copy = copy(obj)
            obj_copy = finiteStateMachine();
        end
    end    
end