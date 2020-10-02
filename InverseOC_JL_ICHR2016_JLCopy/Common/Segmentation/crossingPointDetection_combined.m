function featureModel = crossingPointDetection_combined(templateName, trainingTime, trainingVelocity, trainingLabel, settings)

sigDofArray = [];

for i = 1:length(trainingTime)
    sigDofArray = [sigDofArray; segmentModiftZvcMethod_preamble(trainingVelocity{i}')];
end

sigDof = mode(sigDofArray(sigDofArray>0)); % ignore all zeros
allSigDof = unique(sigDof);

for ind_templateName = 1:length(templateName)
    % now figure out all the labels, and create a model for each
    currTemplateName = templateName{ind_templateName};
    currTrainingTime = {};
    currTrainingVelocity = {};
    
    for ind_trainingTime = 1:length(trainingTime)
        if strcmpi(currTemplateName, trainingLabel{ind_trainingTime})
            currTrainingTime = [currTrainingTime trainingTime(ind_trainingTime)];
            currTrainingVelocity = [currTrainingVelocity trainingVelocity(ind_trainingTime)];
        end
    end
    
    sigDofArray = [];
    for i = 1:length(currTrainingVelocity)
        sigDofArray = [sigDofArray; segmentModiftZvcMethod_preamble(currTrainingVelocity{i}')];
    end
    
    sigDof = mode(sigDofArray(sigDofArray>0)); % ignore all zeros
    
    switch settings
        case 'twoPeak'
            featureModel{ind_templateName} = ...
                crossingPointDetection_twopeak(currTemplateName, currTrainingTime, currTrainingVelocity, sigDof, allSigDof);
            
        case 'multiPeak'
            featureModel{ind_templateName} = ...
                crossingPointDetection_multipeak(currTemplateName, currTrainingTime, currTrainingVelocity, sigDof, allSigDof);
    end
end



% param.allowNeg = 1;
% normPos = normVector(pos(:, sigDof), param);
% normVelo = [0; diff(normPos)];
end

function featureModel = crossingPointDetection_twopeak(templateName, trainingTime, trainingVelocity, sigDof, allSigDof)
    % Feature extraction 

    globalVeloMultiplierAbs = 0.7; 
    globalVeloMultiplier = 0.8; 
    veloThrshold = 2; % TODO threshold
    thresholdObservationLength = 1;
    lenTrainingTime = length(trainingTime);
    allowDoubleZVCs = 0;
    
%     velo = cell(1, lenTrainingTime);
    crossingTime = cell(1, thresholdObservationLength);
    crossingMtx = cell(1, thresholdObservationLength);
    
    sigDofMtx = allSigDof;
    for i = 1:length(allSigDof)
        if isempty(intersect(allSigDof(i), sigDof))
            sigDofMtx(i) = 0;
        end
    end

    for i = 1:thresholdObservationLength

        veloLocal = [];
        timeLocal = [];
        
        crossingMtxAll = cell(lenTrainingTime, 1);
        
        for j = 1:length(trainingTime)
            baseTime = trainingTime{j}; % - trainingTime{j}(1); % start zero
            baseVelo = veloMulti(trainingVelocity{j}, sigDof);
            crossingMtxAll{j} = [];
            
            % determine global max and min
            globalMaxVeloAbs = max(abs(baseVelo));
            globalMaxVelo = max(baseVelo);
            globalMinVelo = min(baseVelo);

            [crossingStruct, crossingDof] = zvcCheck(baseTime, baseVelo);   
            crossingInd = unique([1 crossingStruct{1}.Index length(baseTime)]); % default structure already creates a ZVC at start
            
            doubleZVC(j) = 0; % flag to note if double-zvcs are noticed
            counter = 1;
            counterPeakZVC = 1;
            crossingMtxAll{j}(counterPeakZVC) = 1;
            timeLocal(j, counterPeakZVC) = baseTime(crossingInd(counter));
            
            % examine the velos between each crossing
            crossPeakVelo = zeros(size(crossingInd))';
            for ind_crossing = 2:length(crossingInd)
                veloToExamine = baseVelo(crossingInd(ind_crossing-1):crossingInd(ind_crossing));
                [crossPeakVelo(ind_crossing), crossPeakInd(ind_crossing)] = max(abs(veloToExamine));
                crossPeakInd(ind_crossing) = crossPeakInd(ind_crossing) + crossingInd(ind_crossing-1);
            end
            
            if 0
                figure;
                plot(baseVelo);
                hold on
                scatter(crossingInd, zeros(size(crossingInd)), 'ro');
                
                ack = crossPeakVelo' * 10000;
            end
            
            [maxPeakVal, maxPeakInd] = max(crossPeakVelo);
            crossPeakVelo(maxPeakInd) = 0;
            startZVCInd1 = crossingInd(maxPeakInd-1);
            maxPeakMagInd = crossPeakInd(maxPeakInd);
            startZVCInd2 = crossingInd(maxPeakInd);
            
            [secondPeakVal, secondPeakInd] = max(crossPeakVelo);
            endZVCInd1 = crossingInd(secondPeakInd-1);
            secondPeakMagInd = crossPeakInd(maxPeakInd);
            endZVCInd2 = crossingInd(secondPeakInd);
            
            crossingTime = [];
            zvcCrossingTime{1} = baseTime([startZVCInd1 startZVCInd2 endZVCInd1 endZVCInd2]);
            veloTime{1} = [];
            
            crossingMtx = [];
            sigDofMtx = [];
            veloTime = [];
            veloMag = [];
            peakMag = [];
            peakCount = [];
            crossingVeloInd(1) = 1;
            doubleZvcFlag = [];
        end
        
%         crossingVeloInd = find(crossingMtx{1} ~= 1);
        
        featureModel.name = templateName;
        featureModel.crossingTime = crossingTime;
        featureModel.zvcTime = zvcCrossingTime;
        featureModel.crossingMtx = crossingMtx;
        featureModel.sigDof = sigDofMtx;
        featureModel.veloTime = veloTime;
        featureModel.veloMag = veloMag;
        featureModel.peakMag = peakMag;
        featureModel.peakCount = peakCount;
        featureModel.startingZVCCount = crossingVeloInd(1) - 1;
        featureModel.endingZVCCount = [];
        featureModel.doubleZvcFlag = [];
    end
end

function featureModel = crossingPointDetection_multipeak(templateName, trainingTime, trainingVelocity, sigDof, allSigDof)
    % Feature extraction 

    globalVeloMultiplierAbs = 0.3; 
    globalVeloMultiplier = 0.8; 
    veloThrshold = 2; % TODO threshold
    thresholdObservationLength = 1;
    lenTrainingTime = length(trainingTime);
    allowDoubleZVCs = 0;
    
%     velo = cell(1, lenTrainingTime);
    crossingTime = cell(1, thresholdObservationLength);
    crossingMtx = cell(1, thresholdObservationLength);
    
    sigDofMtx = allSigDof;
    for i = 1:length(allSigDof)
        if isempty(intersect(allSigDof(i), sigDof))
            sigDofMtx(i) = 0;
        end
    end
    
%     for i = 1:lenTrainingTime
%         time1 = trainingTime{i}';
%         data1 = trainingData{i}';
%         
%         velo{i} = veloCalc(time1, data1, 1, filterParam)';
%     end

zvcCrossingTime = [];

    for i = 1:thresholdObservationLength

        veloLocal = [];
        timeLocal = [];
        
        crossingMtxAll = cell(lenTrainingTime, 1);
        
        for j = 1:length(trainingTime)
            baseTime = trainingTime{j}; % - trainingTime{j}(1); % start zero
            baseVelo = veloMulti(trainingVelocity{j}, sigDof);
            crossingMtxAll{j} = [];
            
            % determine global max and min
            globalMaxVeloAbs = max(abs(baseVelo));
            globalMaxVelo = max(baseVelo);
            globalMinVelo = min(baseVelo);

            [crossingStruct, crossingDof] = zvcCheck(baseTime, baseVelo);   
            crossingInd = unique([crossingStruct{1}.Index]); % default structure already creates a ZVC at start
            
            % if there's no ZVCs near the ending, add it as a ZVC
            if length(baseTime) - crossingInd(end) > 10
                 crossingInd = unique([crossingStruct{1}.Index length(baseTime)]); % default structure already creates a ZVC at start
            end
            
            doubleZVC(j) = 0; % flag to note if double-zvcs are noticed
            counter = 1;
            counterPeakZVC = 1;
            crossingMtxAll{j}(counterPeakZVC) = 1;
            timeLocal(j, counterPeakZVC) = baseTime(crossingInd(counter));
            
            while counter < length(crossingInd)
                ind1 = crossingInd(counter);
                counter = counter + 1;
                ind2 = crossingInd(counter);
                absPeakVelo = max(abs(baseVelo(ind1:ind2)));
                
                % the peaks examined has to be of a certain magnitude
                if (max(baseVelo(ind1:ind2)) > globalMaxVelo*globalVeloMultiplier || ...
                        min(baseVelo(ind1:ind2)) < globalMinVelo*globalVeloMultiplier) && ...
                    absPeakVelo > globalMaxVeloAbs*globalVeloMultiplierAbs  % encourages two peak patterns...
%                 if absPeakVelo > globalMaxVeloAbs*globalVeloMultiplierAbs
                    % reached acceptiable thresholds
                    
                    % peak seeker calculate
                    counterPeakZVC = counterPeakZVC + 1;
                    [crossingMtxAll{j}(counterPeakZVC), timeLocal(j, counterPeakZVC), veloLocal(j, counterPeakZVC)] = ...
                        peakSeeker(baseTime, baseVelo, [ind1 ind2]);
                    
                    % there is a ZVC in between the peaks
                    counterPeakZVC = counterPeakZVC + 1;
                    crossingMtxAll{j}(counterPeakZVC) = 1;
                    timeLocal(j, counterPeakZVC) = baseTime(ind2);
                    veloLocal(j, counterPeakZVC) = 0;
                    zvcCrossingTime = [zvcCrossingTime baseTime(ind1) baseTime(ind2)];
                    
                elseif allowDoubleZVCs && (counter == 2 || counter == length(crossingInd))
                    % keep the item as a ZVC if it's first or last entry
                    counterPeakZVC = counterPeakZVC + 1;
                    crossingMtxAll{j}(counterPeakZVC) = 1;
                    timeLocal(j, counterPeakZVC) = baseTime(ind2);
                    veloLocal(j, counterPeakZVC) = 0;
                    
                    doubleZVC(j) = 1;
                else
                    % else discard the ZVC
                    continue
                end
            end
         
            
%             startTime(j) = 0;
%             endTime(j) = baseTime(end);
        end
        
%         testCrossingMtx = crossingMtxAll{1};
        votingMtx = {};
        for j = 1:length(crossingMtxAll)
            matchFlag = 0;
            frontDoubleZVC = 0;
            rearDoubleZVC = 0;
            
            if doubleZVC(j) == 1 
                % checking the front for double zvc
                if isequal(crossingMtxAll{j}(1:2), [1 1])
                    startingInd = 2;
                    frontDoubleZVC = j;
                else 
                    startingInd = 1;
                end
                
                % checking the back for double zvc
                if isequal(crossingMtxAll{j}(end-1:end), [1 1])
                    endingInd = length(crossingMtxAll{j}) - 1;
                    rearDoubleZVC = j;
                else
                    endingInd = length(crossingMtxAll{j});
                end
                
                % noting the actual feature array, without the zvcs
                crossingMtxToMatch = crossingMtxAll{j}(startingInd:endingInd);
                doubleZVC(j) = j;
            else
                crossingMtxToMatch = crossingMtxAll{j};
            end
                
            for k = 1:size(votingMtx, 1)
                if isequal(crossingMtxToMatch, votingMtx{k, 1})
                    % current sequence is already in the matrix
                    votingMtx{k, 2} = [votingMtx{k, 2} j];
                    votingMtx{k, 3} = [votingMtx{k, 3} doubleZVC(j)]; % third entry denotes if double zvc is detected
                    votingMtx{k, 4} = [votingMtx{k, 4} frontDoubleZVC];
                    votingMtx{k, 5} = [votingMtx{k, 5} rearDoubleZVC];
                    matchFlag = 1;
                end
            end
            
            if ~matchFlag
                % make new entry
                votingMtx = [votingMtx; {crossingMtxToMatch j doubleZVC(j) frontDoubleZVC rearDoubleZVC}];
            end
        end
        
        % decide on actual matrix
        electedTemplate = 0;
        currentTopCountTemplate = 0;
        for k = 1:size(votingMtx, 1)
            if length(trainingTime) < 30 % 30 is arb
                % if there is a low template count, go for majority
                if length(votingMtx{k, 2}) >= (2/4)*length(trainingTime)
                    % simple majority
                    electedTemplate = k;
                    break
                end
            else
                % if there is a large template count, just go for most common
                if length(votingMtx{k, 2}) > currentTopCountTemplate
                    % select highest
                    currentTopCountTemplate = length(votingMtx{k, 2});
                    electedTemplate = k;
                end
            end
        end
        
        if electedTemplate > 0
            % selected a good template
        else
            % inconclusive. get the 5 point one
            for k = 1:size(votingMtx, 1)
                if length(votingMtx{k, 1}) == 5
                    electedTemplate = k;
                    break
                end
            end
        end
       
        kcounter = 0;
for k = electedTemplate %1:size(votingMtx, 1)
    %         % now decide on the actual matrix (TODO threshold on the 3/4)
%         if electedTemplate > 0
%             % a majority agrees on the overall template structure
%             k = electedTemplate;

            % else they're mostly/all single zvcs
            doubleZvcEntries = votingMtx{k, 3}(votingMtx{k, 3} > 0);
            entries = setxor(votingMtx{k, 2}, doubleZvcEntries);
            doubleZvcFlag = 0;

            testCrossingMtx = crossingMtxAll{entries(1)}; % allow for double zvc
%                 testCrossingMtx = votingMtx{k, 1}; % no double zvc

            len = length(testCrossingMtx); % just in case the time/velo local array is larger then expected

            timeAvg = mean(timeLocal(entries, 1:len), 1);
            veloAvg = mean(veloLocal(entries, 1:len), 1);

            crossingTime{i} = timeAvg;
            crossingMtx{i} = testCrossingMtx;
            veloTime{i} = timeAvg(abs(veloAvg) > 0);
            veloMag{i} = veloAvg;
            peakMag{i} = veloAvg(abs(veloAvg) > 0);
            peakCount = length(find(abs(testCrossingMtx) == 2));
%         else
%             error('Inconsistent template data');
%         end
        
        crossingVeloInd = find(crossingMtx{1} ~= 1);
        

%         featureModel.zvcTime{1} = zvcCrossingTime;
        featureModel.name = templateName;
        featureModel.crossingTime = crossingTime;
        featureModel.crossingMtx = crossingMtx;
        featureModel.sigDof = sigDofMtx;
        featureModel.veloTime = veloTime;
        featureModel.veloMag = veloMag;
        featureModel.peakMag = peakMag;
        featureModel.peakCount = peakCount;
        featureModel.startingZVCCount = crossingVeloInd(1) - 1;
        featureModel.endingZVCCount = length(crossingMtx{1}) - crossingVeloInd(end);
        featureModel.doubleZvcFlag = doubleZvcFlag;
        
%         kcounter = kcounter + 1;
%         featureModelOut{kcounter} = kcounter;
end
    end
end

function [crossingMtx, timeLocal, veloLocal] = peakSeeker(time, data, zvcPoints)
    % data, assumes to be full
    % zvcpoints, the two pinchoff points (expecting 2 points)
    
    % look at first chunk - expecting either a positive or negative peak
    [maxVeloFront, maxIndFront] = max(data(zvcPoints(1):zvcPoints(end)));
    [minVeloFront, minIndFront] = min(data(zvcPoints(1):zvcPoints(end)));
            
            % which point is closest to the middle?
%             frontInd = [maxIndFront minIndFront];
%             [closeVal, closeInd] = findClosestValue(mean(zvcPoints), frontInd);
            
%             if find(frontInd == closeVal) == 1;
    if abs(maxVeloFront) > abs(minVeloFront) % keep the dominant
        % max is closest to the middle
        crossingMtx = 2;
        timeLocal = time(zvcPoints(1) + maxIndFront - 1); % should have a -1 here
        veloLocal = maxVeloFront;
    else
        % min is closest to the middle
        crossingMtx = -2;
        timeLocal = time(zvcPoints(1) + minIndFront - 1); % should have a -1 here
        veloLocal = minVeloFront;
    end
end

function [crossingStruct, crossingDof] = zvcCheck(timeData, veloToTest, crossingStruct, timeStep, threshold)
    % Check for zero velocity crossings
    
    W = length(veloToTest);
    
    if ~exist('crossingStruct', 'var') || isempty(crossingStruct)
        dof = size(veloToTest, 1);
        
        for i = 1:dof
            
            crossingStruct{i}.Index = [1]; % TODO check this later. Start each crossing logger with a ZVC
            crossingStruct{i}.Time = [0];
            crossingStruct{i}.Mtx = [1];
            crossingStruct{i}.prevCrossingInd = -W;
        end
    end
    
    if ~exist('timeStep', 'var') || isempty(timeStep)
        timeStep = 0; %size(veloToTest, 2);
    end
    
    if ~exist('threshold', 'var')
        threshold = 0;
    end
    
    % window step for ZVC sweeping
%     zvcWindowStep = W-1; % this would be the full size of the W
    zvcWindowStep = 1;
    
%     distanceFromPreviousZVC = floor(W/2);
    distanceFromPreviousZVC = 10; 
    
    crossingDof = zeros(size(veloToTest, 1), 1);
    iInd = 1:size(veloToTest, 1);
    jInd = 1:zvcWindowStep:size(veloToTest, 2)-1;
    
    for i = iInd % each dof...
%         range = max(veloToTest(i, 2:end-1)) - min(veloToTest(i, 2:end-1));
        range = abs(mean(veloToTest(i, :)));
        crossing = zeros(1, length(jInd));

        for j = jInd % 2 timestep window            
            if range < threshold
                % Low velocity entry (change this to 1 and comment out the
                % bottom stuff in order to switch back to previous mode)
                crossing(j) = 0.5; % TODO change to the entire window?
            end
            
            % actual ZVCs is more important
            if veloToTest(i, j) < 0 && veloToTest(i, j+zvcWindowStep) > 0
                % ZVC: if a positive crossing was made
                crossing(j) = 1;
            elseif veloToTest(i, j) > 0 && veloToTest(i, j+zvcWindowStep) < 0
                % ZVC: if a negative crossing was made
                crossing(j) = 1; % was -1 (direction mattered) before %TODO make a decision
            end

            % insert it into the crossing matrix, if it passes the time treshold
            if crossing(j) ~= 0 && timeStep + j >= (crossingStruct{i}.prevCrossingInd + distanceFromPreviousZVC)
                % positive threshold observed
%                 crossingStruct.Counter(i) = crossingStruct.Counter(i) + 1;
    
                nextZvcInd = timeStep + j + floor(zvcWindowStep/2);
                
                if nextZvcInd == 0
                    nextZvcInd = 1;
                elseif nextZvcInd > length(timeData)
                    nextZvcInd = length(timeData);
                end

                currCrossingInd = [crossingStruct{i}.Index nextZvcInd];
                crossingStruct{i}.Index = sort(unique(currCrossingInd)); % feed into timeData to get crossingTime
                crossingStruct{i}.prevCrossingInd = timeStep + j; % update the 'last crossing' item
                
                crossingDof(i) = 1; % the dof that got a crossing
            end
        end
        
%         if timeStep > 120 && isfield(crossingStruct{i}, 'prevCrossingChar') && ...
%                 mean(crossing) == 0.5 && mean(crossingStruct{i}.prevCrossingChar) == 0.5
%             % if this crossing and the previous crossing are close by and
%             % also is a low velo ZVC, combine them to make a 'full crossing'
%             currCrossingInd = [crossingStruct{i}.Index timeStep+ceil(zvcWindowStep/2)];
%             crossingStruct{i}.Index = sort(unique(currCrossingInd)); % feed into timeData to get crossingTime
%             crossingStruct{i}.prevCrossingInd = timeStep - 2*W; % set to the start of the low velo part
%             
%             crossingDof(i) = 1; % the dof that got a crossing
%             
%             crossing = zeros(1, length(jInd)); % erase previous crossing data
%         end

        currCrossingInd = crossingStruct{i}.Index;
        crossingStruct{i}.Time = timeData(currCrossingInd);
        crossingStruct{i}.Mtx = ones(size(currCrossingInd));
        crossingStruct{i}.prevCrossingChar = crossing;
        crossingStruct{i}.Counter = length(currCrossingInd);
    end 

        % check for entry into low velocity zone
%     crossing = zeros(size(veloToTest, 1), 1);
%     for i = 1:size(veloToTest, 1)
%         range = max(veloToTest(i, 2:end-1)) - min(veloToTest(i, 2:end-1));
%         
%         if range < threshold && abs(veloToTest(1)) > threshold / 2
%             crossing(i) = -1;
%         elseif range < threshold && abs(veloToTest(end)) > threshold / 2;
%             crossing(i) = 1;
%         end
%         
% %         if range < threshold && veloToTest(1) > threshold/2
% %             % entering the low region
% %             crossing(i) = -1;
% %         elseif range < threshold && veloToTest(end) > threshold/2
% %             % exiting the low region
% %             crossing(i) = 1;
% %         end
%     end

    % check actual ZVC
%     crossing = zeros(size(veloToTest, 1), 1);  % no crossing was made
%     for i = 1:size(veloToTest, 1)
%         if veloToTest(i, 1) < 0 && veloToTest(i, end) > 0
%             % if a positive crossing was made
%             crossing(i) = 1;
%         elseif veloToTest(i, 1) > 0 && veloToTest(i, end) < 0
%             % if a negative crossing was made
%             crossing(i) = -1;
%         end
%     end
end

function veloM = veloMulti(data, sigDofs)
    veloM = ones(1, size(data, 2));

    for i = 1:length(sigDofs)
        if sigDofs(i) > 0
            veloM = veloM .* data(sigDofs(i), :);
        end
    end
    
%     figure
%     plot(veloM, 'o');
%     hold on
%     plot(data(sigDofs, :)');
end

