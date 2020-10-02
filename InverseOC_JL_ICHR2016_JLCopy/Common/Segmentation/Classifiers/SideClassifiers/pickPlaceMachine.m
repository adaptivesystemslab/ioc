classdef pickPlaceMachine < AClassifier
% FSM
    properties(SetAccess=protected)
        n = 10; % number of training data it takes to make a new point
        ageLimit = 50; % number of datapoints that would be seen for a given class before purging
        overlapCoeff = 0.3;
        
        trainingDataCache = [];
        trainingLabelCache = [];        
        
        trainingDataOrig = [];
        trainingLabelOrig = [];
        
        outputUnique = [];
        clusterCentreArray = [];
        clusterDistanceArray = [];
        clusterLabelArray = [];
        clusterAgeArray = [];
        
        trainingDataPreviousPrepruned = [];
        trainingLabelPreviousPrepruned = [];
        trainingDistancePreviousPrepruned = [];
        
        dimReduct = [];
        
        training = [];
        group = [];
        
        fsm = [];
        
        lastTrainedTimeStamp = [];
    end
    
    methods
        function obj = pickPlaceMachine(incParam)
            obj.n = incParam.datapointsInACluster;
            obj.ageLimit = incParam.ageThresholdForOldDatapoints;
            obj.overlapCoeff = incParam.overlapInDistance;
        end
            
        function init(obj,varargin)
        % Initialization function for FSM

        end
        
        function [clusterCentre, clusterDist, minDist, minInd] = calcKmeans(obj, data, k)
            % calculate kmeans
            [IDX, clusterCentre, SUMD, D] = kmeans(data, k);
            
            dist = sqrt(D); % kmeans distance is reported to be the sq euclidean
            
%             for i = 1:k
%                 distVec(:, i) = normVector(repmat(clusterCentre(i, :), size(data, 1), 1) - data);
%             end
            
            % determine distance
            for i = 1:k
                clusterDist(i) = max(dist(IDX == i, i));
            end
            
            % keep the points with the shortest distance
            for i = 1:k
                [minDist(i), minInd(i)] = min(dist(IDX == i, i));
            end
        end
        
        function updateDimReductClass(obj, dimReduct)
            % save the dim reduction into the class
            obj.dimReduct = dimReduct;
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
                
                if C1 > 0.95
                    % if the min is more than 85% of the way away, ie we're uncertain
                    labelToApply = -1;
%                     probToApply = 1;
%                 elseif (obj.clusterLabelArray(I1) ~= obj.clusterLabelArray(I2) ...
%                         && abs(C2 - C1) < 0.10) 
%                     % two different labels, and are close
%                     labelToApply = -1;
%                     probToApply = 1;
                else
                    labelToApply = obj.clusterLabelArray(I1);
%                     probToApply = 1;
                end
                
                labelArrayOut(i) = labelToApply;
                labelArrayProb(i) = C1;
                
                % increment all the datapoints by age
                obj.clusterAgeArray = obj.clusterAgeArray + 1;
                
                % update the age of the match we got
                obj.clusterAgeArray(I1) = 1;
            end
        end
        
        function error = retrain(obj,input,output)
            % combine the original input with the new output
            
            newInput = [input; obj.trainingDataCache];
            newOutput = [output; obj.trainingLabelCache];
            
            train(obj, newInput, newOutput);
        end
        
        function pruneOldDatapoints(obj)
            % remove all the centroids that are over a certain age
            toKeep = obj.clusterAgeArray < obj.ageLimit;
            
            toKeepInd = find(toKeep);
            
            % if the pruning would remove a given class, keep the youngest
            % points in that class
            for i = 1:length(obj.outputUnique)
                toKeepTest = unique(obj.clusterLabelArray(toKeepInd, :));
                
                if sum(toKeepTest == obj.outputUnique(i)) == 0
                    % the class disappears... 
                    index = 1:length(obj.clusterLabelArray);
                    pointsClass = obj.clusterLabelArray == obj.outputUnique(i);

                    localClusterAgeArray = obj.clusterAgeArray(pointsClass, :);
                    
                    [minVal, minInd] = min(localClusterAgeArray);
                    
                    index = index(pointsClass);
                    actualMinInd = index(minInd);
                    
                    toKeepInd = [toKeepInd; actualMinInd];
                end
            end
            
            
            obj.clusterCentreArray = obj.clusterCentreArray(toKeepInd, :);
            obj.clusterDistanceArray = obj.clusterDistanceArray(toKeepInd, :);
            obj.clusterLabelArray = obj.clusterLabelArray(toKeepInd, :);
            obj.clusterAgeArray = obj.clusterAgeArray(toKeepInd, :);
        end
        
        function resetTrainingData(obj)
%             train(obj, obj.trainingDataOrig, obj.trainingLabelOrig);
            
            % this training will be added to the next retraining
            obj.trainingDataCache = obj.trainingDataOrig;
            obj.trainingLabelCache = obj.trainingLabelOrig;
            
            obj.trainingDataPreviousPrepruned = [];
            obj.trainingLabelPreviousPrepruned = [];
            obj.trainingDistancePreviousPrepruned = [];
        end
        
        function error = train(obj,input,output)
            % add input/output data to the current array
            
            if isempty(obj.lastTrainedTimeStamp)
                % first time training. copy out the original data in case
                % we need to reset the classifier
                obj.trainingDataOrig = input;
                obj.trainingLabelOrig = output;
                obj.outputUnique = unique(output);
            end
            
            % now perform the dim reduction, since the DR may change from
            % instance to instance
            inputDR = obj.dimReduct.apply(input);
            
            % pull out each of the clusters
            diffOutput = find(diff(output));
            pairsStart = [1; diffOutput+1];
            pairsEnd   = [diffOutput(1:end); length(output)];
            pairs = [pairsStart pairsEnd];
            
            counter = 0;
            
            dataToKeep = [];
            outputToKeep = [];
            distanceToKeep = [];
            
            % cluster the NEW data
            for ind_input = 1:length(pairsStart)
                % every n points or label switch, make a new centroid
                for ind_training = pairsStart(ind_input):obj.n:pairsEnd(ind_input)
                    
                    startInd = ind_training;
                    endInd = ind_training+obj.n;
                    
                    if endInd > pairsEnd(ind_input)
                        endInd = pairsEnd(ind_input);
                    end
                    
                    if endInd - startInd > 1
                        % calculate mean and distance
                        counter = counter + 1;
                        
                        dataUse = input(startInd:endInd, :);
                        dataUseDR{counter} = inputDR(startInd:endInd, :);
                        outputToTrain{counter} = unique(output(startInd:endInd, :));
                        
                        [clusterCentre{counter}, clusterDistance{counter}, minDist, minInd] ...
                            = calcKmeans(obj, dataUseDR{counter}, 1);
                        
                        % data to keep
                        dataToKeep = [dataToKeep; dataUse(minInd, :)];
                        outputToKeep = [outputToKeep; outputToTrain{counter}];
                        distanceToKeep = [distanceToKeep; clusterDistance{counter}];
                        
                        %                 elseif numberOfOutput == 1
                        %                     clusterCentre{i} = dataUse{i};
                        %                     clusterDistance{i} = 1;
                        
                        obj.clusterCentreArray = [obj.clusterCentreArray; clusterCentre{counter}];
                        obj.clusterDistanceArray = [obj.clusterDistanceArray; clusterDistance{counter}'];
                        obj.clusterLabelArray = [obj.clusterLabelArray; outputToTrain{counter}];
                        obj.clusterAgeArray = [obj.clusterAgeArray; 1];
                    else
%                         % if there's only one input, we can't normalize for
%                         % distance, so we'll
%                         clusterCentre{counter} = [];
%                         clusterDistance{counter} = [];
                    end
                end
            end
            
            % include any OLD data
            if ~isempty(obj.trainingDataPreviousPrepruned)
                % now perform the dim reduction, since the DR may change from
                % instance to instance
                inputDR = obj.dimReduct.apply(obj.trainingDataPreviousPrepruned);
                outputDR = obj.trainingLabelPreviousPrepruned;
                distanceDR = obj.trainingDistancePreviousPrepruned;
                
                obj.clusterCentreArray = [obj.clusterCentreArray; inputDR];
                obj.clusterDistanceArray = [obj.clusterDistanceArray; distanceDR];
                obj.clusterLabelArray = [obj.clusterLabelArray; outputDR];
                obj.clusterAgeArray = [obj.clusterAgeArray; ones(size(distanceDR))];
                
                dataToKeep = [dataToKeep; obj.trainingDataPreviousPrepruned];
                outputToKeep = [outputToKeep; obj.trainingLabelPreviousPrepruned];
                distanceToKeep = [distanceToKeep; obj.trainingDistancePreviousPrepruned];
            end
            
            % save all the pre-pruned clustered data
            obj.trainingDataPreviousPrepruned = dataToKeep;
            obj.trainingLabelPreviousPrepruned = outputToKeep;
            obj.trainingDistancePreviousPrepruned = distanceToKeep;
            
            obj.trainingDataCache = []; % clear the cache
            obj.trainingLabelCache = [];
            
            % now prune the data
            pruneOverlappingCentroids(obj);
            
            %             % cluster
            if 0
                figure;
                color = {'r', 'g'};
                for i = 1:length(dataUse)
                    hold on
                    plot(dataUse{i}(:, 1), dataUse{i}(:, 2), [color{mod(1, outputToTrain{i})+1} '.']);
                end
                
                for i = 1:size(obj.clusterCentreArray, 1)
                    hold on
                    plot(obj.clusterCentreArray(i, 1), obj.clusterCentreArray(i, 2), [color{mod(1, obj.clusterLabelArray(i))+1} 'x'], 'markersize', 15);
                end

                xlim([-10 10]); ylim([-8 8]);
            end

            obj.lastTrainedTimeStamp = datestr(now, 'yyyy-mm-dd HH-MM-SS');
        end
        
        function pruneOverlappingCentroids(obj)
            % make sure there's no overlap in centroids. prune out the
            % ones that do
            loop = 1;
            overallCounter = 0;
            loopCountMax = length(obj.clusterLabelArray);
            
            while loop
                restartLoop = 0;
                overallCounter = overallCounter + 1;
                currLabels = obj.clusterLabelArray;
                
                distArray = zeros(length(obj.clusterLabelArray));
                toKeep = ones(size(obj.clusterLabelArray));
                
                counter = 0;
                allCombinations = [];
                for ind_centroids1 = 1:length(obj.clusterLabelArray)
                    for ind_centroids2 = ind_centroids1+1:length(obj.clusterLabelArray)
                        counter = counter + 1;
                        allCombinations(counter, 1) = ind_centroids1;
                        allCombinations(counter, 2) = ind_centroids2;
                    end
                end
                
                for ind_centroids = 1:size(allCombinations, 1)
                    % make sure there's no overlap in centroids. prune out the
                    % ones that do
                    ind_centroids1 =  allCombinations(ind_centroids, 1);
                    ind_centroids2 =  allCombinations(ind_centroids, 2);
                    
                    centroid1 = obj.clusterCentreArray(ind_centroids1, :);
                    centroid2 = obj.clusterCentreArray(ind_centroids2, :);
                    range1 = obj.clusterDistanceArray(ind_centroids1);
                    range2 = obj.clusterDistanceArray(ind_centroids2);
                    label1 = obj.clusterLabelArray(ind_centroids1);
                    label2 = obj.clusterLabelArray(ind_centroids2);
                    
                    currDist = normVector(centroid1 - centroid2);
                    totalRange = (range1 + range2) * obj.overlapCoeff;
                    distArray(ind_centroids1, ind_centroids2) = currDist;
                    
                    if currDist > totalRange || label1 == label2
                        % far enough, or have same labels
                        
                    else
                        % it's the only point in a given class. if that's
                        % the case, don't remove it
                        if sum(label1 == obj.clusterLabelArray) > 1
                            toKeep(ind_centroids1) = 0;
                        end
                        
                        if sum(label2 == obj.clusterLabelArray) > 1
                            toKeep(ind_centroids2) = 0;
                        end
                        
                        toKeepInd = find(toKeep);
                        obj.clusterCentreArray = obj.clusterCentreArray(toKeepInd, :);
                        obj.clusterDistanceArray = obj.clusterDistanceArray(toKeepInd, :);
                        obj.clusterLabelArray = obj.clusterLabelArray(toKeepInd, :);
                        obj.clusterAgeArray = obj.clusterAgeArray(toKeepInd, :);
                        restartLoop = 1;   % reset the process since the clusters changed
                    end
                                      
                    if restartLoop 
                        break
                    end
                end
                
                if isequal(currLabels, obj.clusterLabelArray) || ...
                    overallCounter > loopCountMax
                    % exit when
                    loop = 0;
                end
            end
        end
        
        function [out, prob] = classify(obj,input)
            inputDR = obj.dimReduct.apply(input);
            [out, prob] = labelOnDistance(obj, inputDR);
            
            pruneOldDatapoints(obj);
        end
        
        function obj_copy = copy(obj)
            obj_copy = finiteStateMachine();
        end
    end    
end