classdef clusterMachine < AClassifier
% FSM
    properties(SetAccess=protected)
        k = [0 0];
        dataFormatName = []; % which data to use?
        dataFormatDof = [];
        dataFormatQuad = [];
        dataFormatClusterParam = [];
        
        distTol = 0.8;
        
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
        function obj = clusterMachine(incParam)
            obj.dataFormatName = incParam.dataFormatName;
            obj.dataFormatDof = incParam.dataFormatDof;
            obj.dataFormatQuad = incParam.dataFormatQuad;
            obj.dataFormatClusterParam = incParam.clusterCountPerClass;
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
        
        function [labelArrayOut, labelArrayDist] = labelOnDistance(obj, fullDataInput)
            labelArray  = zeros(size(fullDataInput, 1), size(obj.clusterCentreArray, 1));
            labelArrayOut = zeros(size(fullDataInput, 1), 1);
            labelArrayDist = zeros(size(fullDataInput, 1), 1);
            
            % calculate the distance to each point
            for i_ind = 1:size(obj.clusterCentreArray, 1)
                currCentre = obj.clusterCentreArray(i_ind, :);
                currCovar = obj.clusterCovarArray{i_ind};
                
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

                labelArray(:, i_ind) = dataInputDist;
            end
            
            for i_ind = 1:size(labelArray)
%                 [labelDist, labelInd] = min(labelArray(i, :));

                % pull out the min for all classes
                currClassDist = zeros(size(obj.outputUnique));
                for j_ind = 1:length(obj.outputUnique)                    
                    currClassDist(j_ind) = min(labelArray(i_ind, obj.clusterLabelArray == obj.outputUnique(j_ind)));
                end
                
                % distance check
                [labelDist, labelInd] = min(currClassDist);
                distCheckArray = (currClassDist - min(currClassDist))/min(currClassDist);
                distCheckArray = distCheckArray(distCheckArray > 0); % remove the min entry (would be 0 after this)
                
                labelArrayOut(i_ind) = obj.outputUnique(labelInd);
                labelArrayDist(i_ind) = min(distCheckArray); % this allows us to reject ambig points later
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
        
        function error = retrain(obj,dataStructInput, percentageToKeep, distanceToKeep)
            % combine the original input with the new output
            
            if distanceToKeep < 0
                distanceToKeep = obj.distTol; % use the tuned version if not specified                
            end
            
            if ~isempty(dataStructInput)
                % maintain the same number of data points
                numberOfDataToKeep = floor(length(obj.trainingLabel)) - length(dataStructInput.output) - 1;
                
                if numberOfDataToKeep < 0
                    dataToKeep = [];
                else
                    dataToKeep = randperm(length(obj.trainingLabel), numberOfDataToKeep);
                end
                
            else
                numberOfDataToKeep = floor(length(obj.trainingLabel) * percentageToKeep);
                dataToKeep = randperm(length(obj.trainingLabel), numberOfDataToKeep);
            end
            
            if ~isempty(dataStructInput)
                newInput = [obj.trainingData(dataToKeep, :); dataStructInput.input];
                newOutput = [obj.trainingLabel(dataToKeep, :); dataStructInput.output];
                newOutputGndTruth = [obj.trainingLabelGndTruth(dataToKeep, :); dataStructInput.outputGndTruth];
            else
                % if passed in a blank, ignore the state
                newInput = obj.trainingData(dataToKeep, :);
                newOutput = obj.trainingLabel(dataToKeep, :);
                newOutputGndTruth = obj.trainingLabelGndTruth(dataToKeep, :);
            end
            
            fprintf(' (Cluster SD update w/ %u/%u retained dp) \n', length(dataToKeep), length(newOutput));
            
            dataStruct.input = newInput; % the original data
            dataStruct.output = newOutput;
            dataStruct.outputGndTruth = newOutputGndTruth;
            
            train(obj, dataStruct, distanceToKeep);
        end
        
        function resetTrainingData(obj)
%             train(obj, obj.trainingDataOrig, obj.trainingLabelOrig);
            
            obj.trainingData          = obj.trainingDataOrig;
            obj.trainingLabel         = obj.trainingLabelOrig;
            obj.trainingLabelGndTruth = obj.trainingLabelGndTruthOrig;
            
            retrain(obj, [], 1, -1);
        end
        
        function resetStateMachine(obj)
            obj.clusterCentreArray = [];
            obj.clusterDistanceArray = [];
            obj.clusterLabelArray = [];
            obj.clusterCovarArray = [];
        end
        
        function inputDR = applyDataFormat(obj, input, output)
            if ~exist('output', 'var')
                output = [];
            end
            
           switch obj.dataFormatName
                case 'PC'
                    % now perform the dim reduction, since the DR may change from
                    % instance to instance
                    inputDR = obj.dimReduct.apply(input);
                    
                    if obj.dataFormatDof > 0
                        % restrict the DOFs
                        inputDR = inputDR(:, 1:obj.dataFormatDof);
                    end
                    
%                     switch obj.dataFormatQuad
%                         case 1
%                             inputDR(:, 1) =  inputDR(:, 1) - mean(inputDR(:, 1));
%                         otherwise
%                             
%                     end
                    
               case 'PP'
                   dofs = 1:size(input, 2)/2;
                   
                   % PP always only have 2 vectors (q and dq)
                   inputDR(:, 1) = normVector(input(:, dofs));
                   inputDR(:, 2) = normVector(input(:, dofs+max(dofs)));
                   
                    % after the norm, all the dp would be positive, so
                    % removing the mean from dim2 would just add a negative
                    % offset. if we're not norming the dofs together, then
                    % we may need to mean-normalize the other axis as well
                   inputDR(:, 1) =  inputDR(:, 1) - mean(inputDR(:, 1));
                   
%                    switch obj.dataFormatQuad
%                        case 1
%                                 
%                        otherwise
%                            inputDR(:, 1) =  inputDR(:, 1) - mean(inputDR(:, 1));
%                            inputDR(:, 2) =  inputDR(:, 2) - mean(inputDR(:, 2));
%                    end
                   
               case 'CR'
                   % map to a circle                   
                   inputPCA = obj.dimReduct.apply(input);
                   inputDR = mapToCircle(inputPCA, output ,0); % xp = Ppca' * (x-mu) * Pc;
           end
           
           switch obj.dataFormatQuad
               case 1
                   % we only want the data to occupy one quad
                   inputDR(:, 1) = abs(inputDR(:, 1));
                   inputDR(:, 2) = abs(inputDR(:, 2));
                    
               case 4
                   % leave
           end
           
           % if the DOF is 2 or less, we will use the joint angle as the
           % threshold
           if size(inputDR, 2) <= 2
                [polarTheta,polarRho] = cart2pol(inputDR(:, 1),inputDR(:, 2));
                inputDR = polarTheta; 
           end
           
           if 0
               h1 = figure('position', [       1384         147        1240         913]);
               
               h1(1) = subplot(2, 1, 1);
               
               
%                h1 = figure;
               t = 1:length(output);
               plot(inputDR(output == 1, 1), inputDR(output == 1, 2), 'r.'); hold on
               plot(inputDR(output == 0, 1), inputDR(output == 0, 2), 'b.');
               
%                plot(inputDR(:, 1), inputDR(:, 2), 'r.'); hold on
               
               
               h1(2) = subplot(2, 1, 2);
               polar(polarTheta(output == 1), polarRho(output == 1), 'r.'); hold on
               polar(polarTheta(output == 0), polarRho(output == 0), 'b.')
               
%                 polar(polarTheta, polarRho, 'r.');
           end
           
           if 0
%                h1(1) = subplot(2, 1, 1);
h1 = figure; %   plot(t, inputDR(:, 1), 'b.');
               t = 1:size(inputDR, 1);
               plot(t(output == 1), inputDR(output == 1, 1), 'r.'); hold on
               plot(t(output == 0), inputDR(output == 0, 1), 'b.');
               
               h1(2) = subplot(2, 1, 2);
               plot(t(output == 1), polarTheta(output == 1), 'r.'); hold on
               plot(t(output == 0), polarTheta(output == 0), 'b.');
               linkaxes(h1, 'x');
           end
        end
        
        function trainingFDA(obj, inputDR, output)
            % assuming 2 class, 1 k right now
            blah = FDATransform(-1);
            blah.train(inputDR,output);
            inputDRTrans = blah.apply(inputDR);
            
            figure;
            t = 1:length(output); 
            plot(t(output == 0), inputDRTrans(output == 0), 'r.'); hold on
            plot(t(output == 1), inputDRTrans(output == 1), 'g.'); 
            
%              for i = 1:length(obj.outputUnique)
%                  
%              end
        end
        
        function [clusterCentre, clusterDistance, clusterCovar] = trainingClustering(obj, inputDR, output, k)
            resetStateMachine(obj);
            
            for i = 1:length(obj.outputUnique)
                outputSelection = output == obj.outputUnique(i);
                numberOfOutput = length(find(outputSelection));
                dataUse{i} = inputDR(outputSelection, :);
                
                if ~isempty(k) && numberOfOutput > 1
                    % calculate mean and distance
                    [clusterCentre{i}, clusterDistance{i}, clusterCovar{i}] = calcKmeans(obj, dataUse{i}, k(i));
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
        end
        
        function error = train(obj,dataStruct,distanceToKeep)
            if ~exist('distanceToKeep', 'var')
                % distance from core to keep
                distanceToKeep = obj.distTol;
            end
            
            % tuning variables
            switch obj.dataFormatQuad
                case 1
                    kToTry = 1:8;
                    distThresToTry = 0.1:0.1:3.0;
                    
                case 4
                    kToTry = 2:8;
                    distThresToTry = 0.1:0.1:3.0;
            end

            
            input = dataStruct.input;
            output = dataStruct.output;
            outputGndTruth = dataStruct.outputGndTruth;
            
            inputDR = applyDataFormat(obj, input, output);
                  
            % pull out the data for the two classes
            obj.outputUnique = unique(output);
            
            if isempty(obj.lastTrainedTimeStamp)
                % first time training. copy out the original data in case
                % we need to reset the classifier. make sure to save it in
                % the original space and not DR
                obj.trainingDataOrig          = input;
                obj.trainingLabelOrig         = output;
                obj.trainingLabelGndTruthOrig = outputGndTruth;
                
                % also, we need to determine the k count if left ambig
                if obj.dataFormatClusterParam == -1                                        
                    kArray = [];
                    for i_ind = 1:length(kToTry)
                        for j_ind = 1:length(kToTry)
                            kArray(end+1, 1) = kToTry(i_ind);
                            kArray(end, 2) = kToTry(j_ind);
                        end                        
                    end
                    
                    for i_ind = 1:size(kArray, 1)
                        % calculate the cluster mean and covar 
                        [clusterCentre, clusterDistance, clusterCovar] = trainingClustering(obj, inputDR, output, kArray(i_ind, :));
                        [labelArrayOut, labelArrayDist] = classify(obj, input);
                        
                        % penalize for model complexity
                        penModelComplexity(i_ind) = sum(kArray(i_ind, :));
                        
                        % check the distance between the clusters of one
                        % class to another. penalize if models overlap 
                        penModelOverlap(i_ind) = 0;
                        for j_ind = 1:length(clusterCentre)
                            targetCluster = j_ind+1;
                            if targetCluster > length(clusterCentre)
                                targetCluster = j_ind-1;
                            end
                            
                            dataToTest = clusterCentre{targetCluster};
                            for k_ind = 1:size(clusterCentre{j_ind}, 1)
                                dataInputCentred = [repmat(clusterCentre{j_ind}(k_ind, :), size(dataToTest, 1), 1) - dataToTest]';
                                
                                dataInputDist2 = dataInputCentred' * inv(clusterCovar{j_ind}{k_ind}) * dataInputCentred;
                                dataInputDist = diag(dataInputDist2);
                                
                                minDistDiag = diag(clusterCovar{j_ind}{k_ind});
                                
                                if length(minDistDiag) > 2
                                    % care more about the first 2 dims...
                                    minDistDiagUse = mean(minDistDiag(1:2));
                                else
                                    minDistDiagUse = mean(minDistDiag(1:1));
                                end

%                                  minDistDiagUse = mean(minDistDiag);
                                
                                 penModelOverlap(i_ind) = penModelOverlap(i_ind) + sum(dataInputDist < minDistDiagUse);
                            end
                        end
                        
                        kErr(i_ind) = sum(labelArrayOut' == output)/length(output) ...
                            - 2*penModelComplexity(i_ind)/100 - 2*penModelOverlap(i_ind)/100;
                        
%                          plotData(obj, labelArrayOut, labelArrayDist, inputDR, outputGndTruth, clusterCentre, clusterCovar);
                    end
                    
                    % find the best kErr
                    [kVal, kInd] = max(kErr);
                    obj.k = kArray(kInd, :);                    
                else
                    obj.k = [obj.dataFormatClusterParam obj.dataFormatClusterParam];
                end
            end
            
            [clusterCentre, clusterDistance, clusterCovar] = trainingClustering(obj, inputDR, output, obj.k);
            [labelArrayOut, labelArrayDist] = classify(obj, input);
                          
            if isempty(obj.lastTrainedTimeStamp)
                % now determine the best dist threshold to throw away data
                for i_ind = 1:length(distThresToTry)
                    labIndDistGood = labelArrayDist > distThresToTry(i_ind);
                    
                    dErr(i_ind) = sum(labelArrayOut(labIndDistGood)' == output(labIndDistGood))/length(output(labIndDistGood));
                    penDiscard0(i_ind) = 1 - sum(labIndDistGood(output == 0))/length(output(output == 0)); % penalize based on the number of data discarded
                    penDiscard1(i_ind) = 1 - sum(labIndDistGood(output == 1))/length(output(output == 1)); % penalize based on the number of data discarded
                    
                    % if either error discards more than 70% of the points,
                    % ramp up the error
                    if penDiscard0(i_ind) > 0.7 || penDiscard1(i_ind) > 0.7
                        penDiscard0(i_ind) = 1;
                        penDiscard1(i_ind) = 1;
                    end
                end
                
                [dVal, dInd] = max(dErr - 0.5*penDiscard0  - 0.5*penDiscard1);
                obj.distTol = distThresToTry(dInd);
            end
            
            blahplaceholder = 1;
            
            if 0
                plotData(obj, labelArrayOut, labelArrayDist, inputDR, outputGndTruth, clusterCentre, clusterCovar);
            end
            
            % prune the dataset in prep for future runs
            % if the labelArrayProb says it's within n% of a cluster core,
            % let's keep it
            passOn = 1;
            while passOn
                % keep points that are not ambig (boost the amount of points)
                probToKeep = labelArrayDist >= distanceToKeep; 
                outTest = output(probToKeep);
                outTest0 = outTest == 0;
                outTest1 = outTest == 1;
                
                threshold = 500;
                if length(labelArrayDist) < 3000
                    threshold = 50;
                end
                
                if sum(outTest0) > threshold && sum(outTest1) > threshold
                    % we're good
                    passOn = 0;
                else
                    % boost the distance until we have a good probToCount
                    distanceToKeep = distanceToKeep*1.05;
                    fprintf('Not enough points kept. Boosting to %u\n', distanceToKeep);
                end
            end
            
            obj.trainingData = input(probToKeep, :); % save the data in the ORIG space, not DR space
            obj.trainingLabel = output(probToKeep);
            obj.trainingLabelGndTruth = outputGndTruth(probToKeep);
            
            obj.lastTrainedTimeStamp = datestr(now, 'yyyy-mm-dd HH-MM-SS');
        end
        
        function plotData(obj, labelArrayOut, labelArrayDist, inputDR, outputGndTruth, figureName)
            if ~exist('figureName', 'var')
                figureName = [obj.dataFormatName ' - dof' num2str(obj.dataFormatDof) ' - quad' num2str(obj.dataFormatQuad) ...
                ' - cluster' num2str(obj.k(1)) ',' num2str(obj.k(2))];
            end
            
            labelArrayOut = labelArrayOut';
            h = figure('position', [  1375         303         560*2         420]);
            color = {'r', 'g', 'b'};
            t = [1:length(labelArrayOut)]';
            
            for i = 1:length(obj.outputUnique)
                clusterCentre{i} = obj.clusterCentreArray(obj.clusterLabelArray == obj.outputUnique(i), :);
                clusterCovar{i} = obj.clusterCovarArray(obj.clusterLabelArray == obj.outputUnique(i), :);
            end
            
            for i = 1:length(obj.outputUnique)
                subplot(1, 2, 1);
                hold on
                lab = obj.outputUnique(i);
                labInd1 = find(labelArrayOut == lab);
                labIndDistTooFar = labelArrayDist < obj.distTol;
                if size(clusterCentre{i}, 2) > 1
                    plot(inputDR(labInd1, 1), inputDR(labInd1, 2), [color{i} '.'], 'DisplayName', 'ClusterAssignment');
                    plot(inputDR(labIndDistTooFar, 1), inputDR(labIndDistTooFar, 2), ['ko'], 'DisplayName', 'AmbigData');
                else
                    plot(t(labInd1), inputDR(labInd1), [color{i} '.'], 'DisplayName', 'ClusterAssignment');
                    plot(t(labIndDistTooFar), inputDR(labIndDistTooFar), ['ko'], 'DisplayName', 'AmbigData');
                end
                
                title('Cluster');
                
                subplot(1, 2, 2);
                hold on
                labInd2 = find(outputGndTruth == lab);
                
                if size(clusterCentre{i}, 2) > 1
                    plot(inputDR(labInd2, 1), inputDR(labInd2, 2), [color{i} '.'], 'DisplayName', 'GroundTruth');
                else
                    plot(t(labInd2), inputDR(labInd2), [color{i} '.'], 'DisplayName', 'GroundTruth');
                end
                title('Ground');
            end
            
            title(figureName);
            
            %                 legend('show')
            
            for i = 1:length(obj.outputUnique)
                subplot(1, 2, 1);  hold on
                if size(clusterCentre{i}, 2) > 1
                    plot(clusterCentre{i}(:, 1), clusterCentre{i}(:, 2), [color{i} 'x'], 'markersize', 20);
                    for j = 1:length(clusterCovar{i})
                        error_ellipse(clusterCovar{i}{j}(1:2, 1:2), clusterCentre{i}(j, 1:2));
                    end
                else
                    for j = 1:length(clusterCentre{i})
                        plot(t, ones(size(t))*clusterCentre{i}(j), [color{i} 'x']);
                        plot(t, ones(size(t))*clusterCentre{i}(j) + clusterCovar{i}{j}, [color{i}]);
                        plot(t, ones(size(t))*clusterCentre{i}(j) - clusterCovar{i}{j}, [color{i}]);
                    end
                end
                
                subplot(1, 2, 2);  hold on
                if size(clusterCentre{i}, 2) > 1
                    plot(clusterCentre{i}(:, 1), clusterCentre{i}(:, 2), [color{i} 'x'], 'markersize', 20);
                    for j = 1:length(clusterCovar{i})
                        error_ellipse(clusterCovar{i}{j}(1:2, 1:2), clusterCentre{i}(j, 1:2));
                    end
                else
                    for j = 1:length(clusterCentre{i})
                        plot(t, ones(size(t))*clusterCentre{i}(j), [color{i} 'x']);
                        plot(t, ones(size(t))*clusterCentre{i}(j) + clusterCovar{i}{j}, [color{i}]);
                        plot(t, ones(size(t))*clusterCentre{i}(j) - clusterCovar{i}{j}, [color{i}]);
                    end
                end
            end
                
%                 if size(clusterCentre{i}, 2) > 1
%                     subplot(1, 2, 1); xlim([-6 6]); ylim([-6 6]);
%                     subplot(1, 2, 2); xlim([-6 6]); ylim([-6 6]);
%                 end
                
                saveas(h, ['C:\Documents\MATLABResults\imgDump\' datestr(now, 'yyyy-mm-dd-HH-MM-SS') '_' figureName '.jpg']);
                saveas(h, ['C:\Documents\MATLABResults\imgDump\' datestr(now, 'yyyy-mm-dd-HH-MM-SS') '_' figureName '.fig']);
                close(h);
        end
        
        function [out, dist] = classify(obj,input)
            jumps = 1000;
            for i = 1:jumps:size(input, 1) % memory overflow issues when the array is too big
                indStart = i;
                indEnd = i+jumps;
                if indEnd > size(input, 1)
                    indEnd = size(input, 1);
                end
                
                inputDR = applyDataFormat(obj, input(indStart:indEnd, :));
                [out(indStart:indEnd), dist(indStart:indEnd)] = labelOnDistance(obj, inputDR);
            end
        end
        
        function classifyAndPlot(obj, input, dataStruct)
            
            [labelArrayOut, labelArrayDist] = classify(obj, input);
            
%             input = dataStruct.input;
%             output = dataStruct.output;
            outputGndTruth = dataStruct.outputGndTruth;
            figureName = dataStruct.name;
            
            inputDR = applyDataFormat(obj, input);
            plotData(obj, labelArrayOut, labelArrayDist, inputDR, outputGndTruth, figureName);
       
        end
        
        function obj_copy = copy(obj)
            obj_copy = finiteStateMachine();
        end
    end    
end