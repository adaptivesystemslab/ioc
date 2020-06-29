function [segResults, observation, fhmmIntermediate, algSegmentResult] = fgSegment(observation, templateStruct, sysparam, fhmmIntermediate)
    % call this function from the Java interface. updates the system with
    % new data array entries for the feature-HMM segmentation

%     global templateStruct observation fhmmIntermediate sysparam
 
    % update the observation array with new data from the Java side
    
   
    [algSegmentResult, observation, fhmmIntermediate] = hmmSegmentFct(templateStruct, observation, fhmmIntermediate, sysparam);
    
    segResults = algSegmentResult.segmentTime';
end
%     observation.currDataLength = 0;
%     observation.currDataChecked = 0;

%% Feature Point Calibrations
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
    
    if ~exist('timeStep', 'var')
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
            if (crossing(j) == 1 ) && timeStep + j >= (crossingStruct{i}.prevCrossingInd + distanceFromPreviousZVC)
                % positive threshold observed % || mean(crossing) == 0.5
%                 crossingStruct.Counter(i) = crossingStruct.Counter(i) + 1;
    
                if crossing(j) == 1
                    nextZvcInd = timeStep + j + floor(zvcWindowStep/2);
                elseif mean(crossing) == 0.5
                    nextZvcInd = timeStep + floor(zvcWindowStep/2);
                end
                
                if nextZvcInd == 0
                    nextZvcInd = 1;
                elseif nextZvcInd > length(timeData)
                    nextZvcInd = length(timeData);
                end

                currCrossingInd = [crossingStruct{i}.Index nextZvcInd];
                crossingStruct{i}.Index = sort(unique(currCrossingInd)); % feed into timeData to get crossingTime
                crossingStruct{i}.prevCrossingInd = timeStep + j; % update the 'last crossing' item
                
                ack = timeData(crossingStruct{i}.Index);
                
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

%% Observation data calibration
function output = filterData(filterVal, filterType, filterParam)
    % expects data in vertical form (each col is a dof)
    
    switch filterType
        case 'noFilter'
            % no filtering
            output = filterVal;
         
        case 'singlePass'
            % Use a singlepass Butterworth
            if exist('filterParam', 'var')
                output = filter_singlepassBW(filterVal, filterParam.bwFreq, filterParam.bwSample, filterParam.bwOrder);
            else
                output = filter_singlepassBW(filterVal);
            end
            
        case 'doublePass'
            % Use a singlepass Butterworth
            if exist('filterParam', 'var')
                output = filter_dualpassBW(filterVal, filterParam.bwFreq, filterParam.bwSample, filterParam.bwOrder);
            else
                output = filter_dualpassBW(filterVal);
            end
            
        case 'movingAverage'
            % moving average filter
            output = zeros(size(filterVal));
            for i = 1:size(filterVal, 2)
                % Use a N-step moving average
                output(:, i) = filter_LPFMovingAverage(filterVal(:, i), filterParam.maWindow);
            end
    end
end

% function [output, zf] = filterData(filterVal, 

%% HMM Algorithms
function [algSegmentResult, observation, fhmmIntermediate] = hmmSegmentFct(templateStruct, observation, fhmmIntermediate, sysparam)
    % Segmentation function
    [algSegmentResult, observation, fhmmIntermediate] = hmmLLMatching_featureWindow_multipeak(templateStruct, observation, fhmmIntermediate, sysparam);
end

function [algSegmentResult, observation, fhmmIntermediate] = hmmLLMatching_featureWindow_multipeak(templateModel, observation, fhmmIntermediate, sysparam)
    % find segmentation points based on LL
    % expect a multipeak system (which can be two peaks, of course)
    
    loglikMtx = [];
    segmentTime = [];
    segmentId = [];
    
	% copying out threshold and setting values to local variables for ease of use
    veloMagThreshold = 2;
    timeThresholdPassFlag = 0;
    
    LLThresholdOverall = templateModel.hmmModel{1}.LLThresholdAbs;
    
    % setting variables
    nSpareZVC = sysparam.segment.NspareZVC; % (was 1 earlier) how much ZVCs will we look ahead for the terminating point?
    thresholdObservationLength = 0.8; % (actually...right now it's checking for magnitude...) TODO should probably make this template specific as well (threshold for merging ZVCs in second sweep)

    lenRecycleTime = sysparam.segment.recycleTime;

    thresTemplateTime = sysparam.feature.lengthMultiplier*templateModel.minTemplateTime;
    thresPeakToPeakTime = sysparam.feature.peaktopeakMultiplier*(templateModel.meanPeakToPeak);
    lookaheadZVCThreshold = sysparam.feature.lookaheadZVC;
    lookaheadStartThreshold = 3; % time 

%     templateLen = 150; %floor(templateModel.meanTemplateLength);
    
    W = sysparam.segment.Wzvc;
    obsTime = observation.obsTime(1:observation.currDataLength);
    obsData = observation.obsData(:, 1:observation.currDataLength);
    
    dof = size(obsData, 1);
    
%     hmmModel = templateModel.hmmModel;
    featureModel = templateModel.featureModel;
    allSigDof    = templateModel.allSigDof;
    lenTemplate  = templateModel.lenTemplate; % TODO totally clean up this section
    lenAllSigDof = length(allSigDof);
        
    for i = 1:lenTemplate
        LLThreshold(i, 1) = templateModel.hmmModel{i}.LLThreshold; % TODO template
    end
    
    % assume starting from rest
    obsVelo = observation.obsVelo(:, 1:observation.currDataLength);
    obsVeloM = ones(lenTemplate, size(obsData, 2));
    
    % create a obsVeloM branch for every template
    for i = 1:lenTemplate
        obsVeloM(i, :) = veloMulti(obsVelo, featureModel{i}.sigDof);
    end
    
    % worker variables
    templateZero = zeros(lenTemplate, 1);
    
    % set up the system state machine
    if observation.currDataChecked == 0
        stateStruct.runHMM = 0; % when 1, start the HMM sequence on the queue
        stateStruct.templateActivated = templateZero;
        stateStruct.waitForNextZVC = templateZero;
        stateStruct.tprevSegment = 1; % marks the end of the previous segment
        
        clearHMMQueue.template = 0;
        clearHMMQueue.state = '';
        clearHMMQueue.pivotPeak = [];
        clearHMMQueue.pivotPeakMag = [];
        clearHMMQueue.edgeEdgeZVCTime = 0;
        clearHMMQueue.endEdgeZVCTime = 0;
        clearHMMQueue.currentZVCLen = 0;
        clearHMMQueue.prePeakZVCInd = [];
        clearHMMQueue.postPeakZVCInd = [];
        
        segmentStruct.segmentId = [];
        segmentStruct.segmentTime = [];
        segmentStruct.segmentMinLL = [];
        segmentStruct.loglikMtx = zeros(length(obsData), 1);
        segmentStruct.segmentCounter = 0;
        
        crossingStruct = cell(1, lenTemplate);
        veloStruct = cell(1, lenTemplate);
        for i = 1:lenTemplate
            crossingStruct{i}.dof = featureModel{i}.sigDof;
            crossingStruct{i}.Index = [];
            crossingStruct{i}.Time = [];
            crossingStruct{i}.Mtx = [];
            crossingStruct{i}.prevCrossingInd = -W;
            crossingStruct{i}.Counter = 0;
            
            veloStruct{i}.dof = featureModel{i}.sigDof;
            veloStruct{i}.Index = [];
            veloStruct{i}.IndexZVCLastChecked = 1;
            veloStruct{i}.IndexVeloLastChecked = 1;
            veloStruct{i}.Time = [];
            veloStruct{i}.Mtx = [];
            veloStruct{i}.Mag = [];
            veloStruct{i}.lastAttenuated = [];
            veloStruct{i}.Counter = 0;
        end
        
        for i = 1:lenTemplate
            stateStruct.hmmQueue(i) = clearHMMQueue;
        end
        
        observation.currDataChecked = 1;
    else
        stateStruct    = fhmmIntermediate.stateStruct;
        clearHMMQueue  = fhmmIntermediate.clearHMMQueue;
        crossingStruct = fhmmIntermediate.crossingStruct;
        veloStruct     = fhmmIntermediate.veloStruct;
        segmentStruct  = fhmmIntermediate.segmentStruct;
    end

    for t = observation.currDataChecked:floor(W/2):size(obsData, 2)-W  
        % decide what to test
        windowedData = obsData(allSigDof, t:t+W);
        windowedTime = obsTime(t:t+W);
        windowedVelo = obsVeloM(:, t:t+W);
        
        if 0
            clf
            subplot(211);
            plot(obsTime, obsData);
            hold on
            plot(windowedTime, windowedData, 'ro');
            
            subplot(212);
            plot(obsTime, obsVeloM);
            hold on
            plot(windowedTime, windowedVelo, 'ro');
        end
        
        % adding to ZVCs
        [crossingStruct, crossingDof] = zvcCheck(obsTime, windowedVelo, crossingStruct, t, sysparam.template.zvcThreshold);
        veloStruct = veloCheck(obsTime, obsVeloM, t, veloStruct, crossingStruct);
        veloStruct = attenuateVelo(veloStruct, windowedTime(1), sysparam);
        timeSinceLastSegment = obsTime(t) - obsTime(stateStruct.tprevSegment);
        
        clearDof = zeros(1, lenAllSigDof);
           
        % jump window time
        if min(windowedTime) > 2.20
            asfsad = 1;
        end
        
        if min(windowedTime) > 50
            asfsad = 1;
        end
        
        if timeSinceLastSegment > thresTemplateTime
            for j = 1:lenTemplate
                passFlag = 0; 
                
                if stateStruct.templateActivated(j) == 0
                    [passFlag, peakInd, veloStruct{j}, pivotPeakTime, pivotPeakMag] = checkVeloStruct(featureModel{j}, crossingStruct{j}, veloStruct{j}, sysparam);
                end
                
                if passFlag == 1
                    veloStruct{j}.IndexVeloPass = peakInd;
                    
                    % add an entry to the HMM queue
                    stateStruct.templateActivated(j) = windowedTime(1);
                    stateStruct.hmmQueue(j).template = j;
                    stateStruct.hmmQueue(j).state = 'waitZVC'; % wait for the ZVC queues to fill up
                    stateStruct.hmmQueue(j).pivotPeak = pivotPeakTime;
                    stateStruct.hmmQueue(j).pivotPeakMag = pivotPeakMag;
                    stateStruct.hmmQueue(j).endEdgeZVCTime = windowedTime(1); % this would be used to determine how long pass the ZVC threshold to wait
                    stateStruct.hmmQueue(j).currentZVCLen = length(crossingStruct{j}.Index);
                end
            end
        end
        
        % check for ZVC after passflags
        for j = 1:lenTemplate
             if stateStruct.templateActivated(j) > 0 ...
                     && strcmpi(stateStruct.hmmQueue(j).state, 'waitZVC')
                 
                 passFlag = 0;
                 if length(crossingStruct{j}.Index) >= nSpareZVC + stateStruct.hmmQueue(j).currentZVCLen
                     % correct number of ZVCs
                     passFlag = 1;
                 end
                 
                 if windowedTime(1) > stateStruct.hmmQueue(j).endEdgeZVCTime + lookaheadZVCThreshold
                     % proper amount of time has elapsed
                     passFlag = 1;
                 end
                 
                 if t > length(obsData)-3*W
                     % at the end of the testing sequence. no other data is
                     % left to window, but we still have unconsumed ZVCs
                     passFlag = 1;
                 end
                 
                 if passFlag
                     stateStruct.hmmQueue(j).state = 'waitHMM';
                     stateStruct.waitForNextZVC(j) = 0; % if we were waiting for an ZVC, clear it
                 end
                 
                 % check for other dofs to see if there are activations. if 
                 % there are large velocities in the other dofs (only if 
                 % those dofs have templates with has that dof as sigdof), 
                 % queue the activetemplate, until all active templates are 
                 % filled
                 if passFlag && t < length(obsData)-3*W && ~isempty(find(featureModel{j}.sigDof == 0))

                     % examine all the dofs labeled zero to see how active
                     % which templates correspond to that?
                     otherTemplates = templateModel.allSigDofExemplar{find(featureModel{j}.sigDof == 0)};
                     
                     for jk = 1:length(otherTemplates)
                         if ~isempty(veloStruct{jk}.Mag) ...
                                 && max(abs(veloStruct{jk}.Mag)) > featureModel{jk}.peakMag{1}(1) * sysparam.feature.velocityMultiplierThreshold ...
                                 && stateStruct.waitForNextZVC(jk) == 0 && stateStruct.templateActivated(jk) == 0
                             % there's a velocity in this template that we
                             % should wait on...
                             stateStruct.waitForNextZVC(jk) = windowedTime(1);
                         end
                     end
                 end
             end
        end

        if sum(stateStruct.templateActivated) > 0
            stateStruct.runHMM = 1;
        end
           
        
        if 0 %min(stateStruct.templateActivated(stateStruct.templateActivated > 0)) < windowedTime(1) + lookaheadStartThreshold
            % if it has been too long since a passflag has been triggered...
            stateStruct.runHMM = 1;
        else
            % otherwise, see if the other templates are we're waiting for
            % have hit the extra ZVCs that we were waiting for
            for j = 1:lenTemplate
                % if all the entries says 'waitHMM' and we've collected all
                % the ZVCs we're looking for, or some time treshold (as
                % determined by the ZVC lookahead threshold) was hit, proceed
                % with the code

                if stateStruct.templateActivated(j) > 0 && strcmpi(stateStruct.hmmQueue(j).state, 'waitZVC')
                    % if any of the queued templates are still waiting for ZVC,
                    % don't proceed with the HMM
                    stateStruct.runHMM = 0;
                end
                
                % if this particular template is waiting for a ZVC...
                if stateStruct.waitForNextZVC(j) > 0
                    if isempty(crossingStruct{j}.Time) || max(crossingStruct{j}.Time) < stateStruct.waitForNextZVC(j) ...
                            || ~(windowedTime(1) > stateStruct.waitForNextZVC(j) + lookaheadZVCThreshold || t > length(obsData)-3*W)
                        % if no new ZVC was added or the threshold time hasn't
                        % been reached, don't proceed with the HMM
                        stateStruct.runHMM = 0;
                    end
                end
            end
        end
        
        % if we're still going to run the HMM...
        if stateStruct.runHMM == 1
            % form each of the queued entry for processing 
            
            stateStruct.firstPeakArray = [];
            stateStruct.secondPeakArray = [];
            
            for j = 1:lenTemplate
                if stateStruct.templateActivated(j) > 0
                    currentFeatureCrossing = featureModel{j}.crossingMtx{1};
                    currentCrossingMtx = unique(crossingStruct{j}.Mtx); % this is '1'
                    currentCrossingTime = crossingStruct{j}.Time;
                    currentCrossingIndex = crossingStruct{j}.Index;
                    pivotPeakTime = stateStruct.hmmQueue(j).pivotPeak;
                    
                    zvcIndexUse = [];
                    zvcIndexUse2 = [];
                    
                    for k = 1:length(pivotPeakTime)
                        currentFirstVelo = pivotPeakTime{k}(1);
                        currentSecondVelo = pivotPeakTime{k}(end);
                        
                        stateStruct.firstPeakArray = [stateStruct.firstPeakArray currentFirstVelo];
                        stateStruct.secondPeakArray = [stateStruct.secondPeakArray currentSecondVelo];
                        
                        [zvcIndex, zvcTime] = zvcTimeReturn_multipeak(obsTime, currentFeatureCrossing, ...
                            currentCrossingMtx, currentCrossingTime, currentCrossingIndex, currentFirstVelo, currentSecondVelo, 'first');
                        
                        [zvcIndex2, zvcTime2] = zvcTimeReturn_multipeak(obsTime, currentFeatureCrossing, ...
                            currentCrossingMtx, currentCrossingTime, currentCrossingIndex, currentFirstVelo, currentSecondVelo, 'last');
                        
%                         % filter out undesiable ZVCs
%                         zvcIndexUse = zvcIndex;
%                         for zvcI = 2:length(zvcIndex)
%                             frontEndData = obsVeloM(j, zvcIndex(zvcI-1):zvcIndex(zvcI));
% %                 if abs(max(frontEndData) - min(frontEndData)) < thresholdObservationLength
%                             if max(abs(frontEndData)) < localThreshold && featureModel{j}.doubleZvcFlag == 0
%                                 % if the range between these two ZVC points are small,
%                                 % then remove the earlier one
%                                 zvcIndexUse(zvcI-1) = 0;
%                             end
%                         end
%                         stateStruct.hmmQueue(j).prePeakZVCInd = zvcIndexUse(zvcIndexUse > 0);
%                         
%                         zvcIndex2Use = zvcIndex2;
%                         for zvcI = 1:length(zvcIndex2)-1
%                             backEndData = obsVeloM(j, zvcIndex2(zvcI):zvcIndex2(zvcI+1));
% %                 if abs(max(backEndData) - min(backEndData)) < thresholdObservationLength
%                             if max(abs(backEndData)) < localThreshold && featureModel{j}.doubleZvcFlag == 0
%                                 % if the range between these two ZVC points are small,
%                                 % then remove the later one
%                                 zvcIndex2Use(zvcI+1) = 0;
%                             end
%                         end
                        
                        % pull the closest ZVCs to the peaks out, up to the
                        % 'allowance' constant
                        
                        % pull the pre-firstpeak ZVCs
                        endingInd = length(zvcIndex) + 1 - featureModel{j}.startingZVCCount;
                        
                        if featureModel{j}.doubleZvcFlag == 1
                            startingInd = endingInd - nSpareZVC * featureModel{j}.startingZVCCount;
                        else
                            startingInd = endingInd;
                        end
                        
                        if startingInd < 1
                            startingInd = 1;
                        end
                        
                        zvcIndexUse = [zvcIndexUse zvcIndex(startingInd:featureModel{j}.startingZVCCount:endingInd)];
                        zvcTimeDebug = zvcTime(startingInd:endingInd);
                        
                        % pull the post-lastpeak ZVCs
                        startingInd = featureModel{j}.endingZVCCount;
                        if featureModel{j}.doubleZvcFlag == 1
                            endingInd = startingInd + nSpareZVC * featureModel{j}.endingZVCCount;
                        else
                            endingInd = startingInd;
                        end
                        
                        if endingInd > length(zvcIndex2)
                            endingInd = length(zvcIndex2);
                        end
                        
                        zvcIndexUse2 = [zvcIndexUse2 zvcIndex2(startingInd:featureModel{j}.endingZVCCount:endingInd)];
                        zvcTime2Debug = zvcTime2(startingInd:endingInd);
                    end

                    stateStruct.hmmQueue(j).prePeakZVCInd = zvcIndexUse(zvcIndexUse > 0);
                    stateStruct.hmmQueue(j).postPeakZVCInd = zvcIndexUse2(zvcIndexUse2 > 0);
                end
            end
        end

		% Start HMM sequence
        if stateStruct.runHMM == 1               
            
            % assmeble the bounds and templates to use
            zvcIndToTest = [];
            hmmTemplateToTest = [];
%             augmentArrayZVC = []; % moved to another line (about 10 lines down)
            minZVCClearTime = inf; % this value is for in case none of the ZVC windowing here are of the appropriate size, recycletime needs to be defined and HMM is skipped here
            peakVeloArrayStuff = [];
            
            for j = 1:lenTemplate
                if stateStruct.templateActivated(j) > 0
                    prePeakLength = length(stateStruct.hmmQueue(j).prePeakZVCInd);
                    postPeakLength = length(stateStruct.hmmQueue(j).postPeakZVCInd);
                    
%                     augmentArrayZVC = zeros(2, prePeakLength*postPeakLength);
                    counter = 0;
                    augmentArrayZVC = [];
                    peakVeloArrayCurrent = [];
                    
                    for j1 = 1:prePeakLength
                        for j2 = 1:postPeakLength
                            if stateStruct.hmmQueue(j).postPeakZVCInd(j2) <= stateStruct.hmmQueue(j).prePeakZVCInd(j1)
                               continue 
                            end
                            
                            prePeakInd = stateStruct.hmmQueue(j).prePeakZVCInd(j1);
                            postPeakInd = stateStruct.hmmQueue(j).postPeakZVCInd(j2);
                            ZVCwindowTime = obsTime(postPeakInd) - obsTime(prePeakInd);
                            templateAvgLength = featureModel{j}.crossingTime{1}(end);
                            
                            if minZVCClearTime > obsTime(prePeakInd);
                                minZVCClearTime = obsTime(prePeakInd);
                            end
                            
                            % check the length, want it within this threshold
                            if ZVCwindowTime > templateAvgLength*(sysparam.feature.windowThresholdMin) ...
                                    && ZVCwindowTime < templateAvgLength*(sysparam.feature.windowThresholdMax)
                            
                                counter = counter + 1;
                                augmentArrayZVC(1, counter) = prePeakInd;
                                augmentArrayZVC(2, counter) = postPeakInd;

                                peakVeloArrayCurrent = [obsVeloM(j, prePeakInd); ...
                                    stateStruct.hmmQueue(j).pivotPeakMag{1}(1)'; ...
                                    obsVeloM(j, postPeakInd)];
                            else
                                if ZVCwindowTime < templateAvgLength*(sysparam.feature.windowThresholdMin) 
                                    fprintf('S');
                                else
                                    fprintf('L');
                                end
                            end
                        end
                    end
                    
                    peakVeloArrayStuff = [peakVeloArrayStuff peakVeloArrayCurrent];
                    zvcIndToTest = [zvcIndToTest augmentArrayZVC];
                    hmmTemplateToTest = [hmmTemplateToTest ones(1, size(augmentArrayZVC, 2))*j];
                end
            end
            
            if ~isempty(hmmTemplateToTest)
                loglikMtxSeg = [];
                loglikMtxDof = [];
                windowSize = [];
                
                for n = 1:length(hmmTemplateToTest)
                    segStart = zvcIndToTest(1, n);
                    segEnd = zvcIndToTest(2, n);
                    windowSize(n) = segEnd - segStart;
                end
                
                templateLen = min(windowSize);
                
                for n = 1:length(hmmTemplateToTest)
                    segStart = zvcIndToTest(1, n);
                    segEnd = zvcIndToTest(2, n);
                    
                    windowedData = obsData(:, segStart:segEnd);
%                 figure; plot(windowedData');
%                 windowSize(n) = 0; % templateLen? segEnd-segStart?
                    windowedDataSpline = resizeObservationData(windowedData, templateLen);
%                 windowedDataSpline = windowedData;
                    
       
        
%         windowToCheck = windowedDataSpline(1:5, :);

%         fprintf('lala\n');
%          templateModel.hmmModel{1}
%         size(windowToCheck)
        
                    % pass in window to check LL
%                     [loglikMtxSeg(n), loglikMtxDof(n)] = hmmLogLikelihood(templateModel, windowToCheck, [], hmmTemplateToTest(n));
                    loglikMtxSeg(n) = 100;
                    loglikMtxDof(n) = hmmTemplateToTest(n);
                end


%     hmmTemplate_time_windowSize_LL = [hmmTemplateToTest; debugTime; windowSize; loglikMtxSeg]
            
            [loglik, loglikInd, templateId] = maxLL(loglikMtxSeg, loglikMtxDof, LLThreshold, LLThresholdOverall, windowSize);
            maxStartInd = zvcIndToTest(1, loglikInd);
            maxEndInd = zvcIndToTest(2, loglikInd);
            maxStartTime = obsTime(maxStartInd);
            maxEndTime = obsTime(maxEndInd);
            
            % pull the other loglikInd that has the same start/end ind
            overlapStartInd = find(zvcIndToTest(1, :) == maxStartInd);
            overlapEndInd = find(zvcIndToTest(2, :) == maxEndInd);
            overlapBoth = intersect(overlapStartInd, overlapEndInd);
            
%             segmentMinLLcurrent = 0;
%             segmentMinLL = 0;

            segmentMinLLvectorTemp = loglikMtxSeg(overlapBoth);
            segmentMinDofMapping = loglikMtxDof(overlapBoth);

            segmentMinLLvector = NaN([lenTemplate 1 1]);
            segmentMinLLvector(segmentMinDofMapping) = segmentMinLLvectorTemp;

%             segmentLLother = setxor(segmentMinLLvector, loglik);  
            segmentLLotherSum = sum(segmentMinLLvector > LLThreshold);
            
            templateMatch = 0;
            
            % break point marker !!!
            % -----------------------------------------------------------------
            if loglik > LLThreshold(templateId) && loglik > LLThresholdOverall
%             if (loglik > LLThreshold(templateId) && ~isempty(intersect(featureModel{templateId}.sigDof, activeSigDof))) || ...
%                     (loglik > LLThreshold(templateId) && segmentLLotherSum < size(loglikMtxSeg, 3))
                
                % 1. is the main LL larger than threshold and...
                % 2. are the other LLs smaller than threshold, allowing for a 'clear winner'?
                % 3. the template that won, is it of the correct dof?
                if sysparam.verbose
                fprintf([num2str(maxStartTime, '%10.2f'), ':' num2str(maxEndTime, '%10.2f'), ' - ']);
                fprintf(['SEGMENT DETECTED - template #', num2str(templateId), ' (', templateModel.motionName{templateId}, ') - LL = ' num2str(loglik) '\n']);
%                 fprintf(['  Peak velo: ...' ...
%                     num2str(peakVeloArrayStuff(1, loglikInd), 3) ' - ' ...
%                     num2str(peakVeloArrayStuff(2, loglikInd), 3) ' - ' ...
%                     num2str(peakVeloArrayStuff(3, loglikInd), 3) ' - ' ...
%                     num2str(peakVeloArrayStuff(4, loglikInd), 3) '\n']);
                
%                 fprintf(['  Peak velo: ...' ...
%                     num2str(peakVeloArrayStuff(1, loglikInd), 3) '\n']);
                end
                
                segmentStruct.segmentCounter = segmentStruct.segmentCounter + 1;
                segmentStruct.segmentTime(1, segmentStruct.segmentCounter) = maxStartInd; %maxStartTime;
                segmentStruct.segmentTime(2, segmentStruct.segmentCounter) = maxEndInd; %maxEndTime;
                segmentStruct.segmentId(segmentStruct.segmentCounter) = templateId;
                
                recycleTime = maxEndTime;
                
%                 if ~isequal(featureModel{templateId}.sigDof, activeSigDof)
%                     fprintf('sigdof mismatch\n');
%                 end
                
                % clear if there is a match
                templateMatch = 1;
                
            else
                % not within threshold range...
                fprintf([num2str(maxStartTime, '%10.2f'), ':' num2str(maxEndTime, '%10.2f'), ' - ']);
                if loglik < LLThreshold(templateId) || loglik < LLThresholdOverall
                    % copied from twopeak
                    if sysparam.verbose
                    fprintf('not over threshold... (is %5.2f, but thres for %s = %5.2f) \n', loglik, templateModel.motionName{templateId}, LLThreshold(templateId));
                    end
                    
                    recycleTime = max(stateStruct.secondPeakArray);
                    
                    % the original one from multipeak
%                     fprintf(['not over threshold for template ' num2str(activeTemplate) '\n']);
%                     recycleTime = min(min([currentFirstVelo currentSecondVelo]));
                    
                elseif segmentLLotherSum > 1
                    if sysparam.verbose
                    fprintf('other LLs are above threshold... \n');
                    end
                    
                    recycleTime = min(stateStruct.firstPeakArray);
%                     recycleTime = min(crossingMtx{activeSigDofInd});
                    
                elseif isempty(intersect(featureModel{templateId}.sigDof, activeSigDof))
                    if sysparam.verbose
                    fprintf('minLL template does not intersect with initiated sigdof \n');
                    end
                    
                    recycleTime = min(stateStruct.firstPeakArray);
                    
                end
            end
            
            stateStruct.tprevSegment = maxEndInd;
            segmentStruct.loglikMtx(maxEndInd) = loglik;
            
            else
                % empty HMMtoTest. None of the ZVC looked at was appropriately sized

                templateMatch = 0;
                recycleTime = minZVCClearTime;
                if sysparam.verbose
                    fprintf('all windows are either too large or too small, skipping one HMM proposal \n');
                end
            end

            % erase the other DOFs that have segment points before the
            % earliest on listed here

                if sysparam.segment.overlapExemplar == 0 && templateMatch
                    % there are no overlapping exemplars, and there is a 
                    % sucessfuly matchclear all the dofs
                    templatesToClear = 1:length(crossingStruct);
                else 
                    % there are overlapping exemplars, clear only the
                    % active ones
                     templatesToClear = 1:length(crossingStruct);
                     % still clear them all, but make the clearing field
                     % smaller
                     recycleTime = min(stateStruct.firstPeakArray);
                     
%                     templatesToClear = templateModel.allSigDofExemplar{find(featureModel{templateId}.sigDof > 0)};
                end
                
                for j = templatesToClear
                    % clean the crossing and velocity structs
                    [crossingStruct{j}, veloStruct{j}] = cleanCrossingStructs(crossingStruct{j}, veloStruct{j}, recycleTime-lenRecycleTime);
                
                    % clean the state machine
                    stateStruct.runHMM = 0;
                    stateStruct.templateActivated(templatesToClear) = 0;
                    stateStruct.waitForNextZVC(templatesToClear) = 0;
                    stateStruct.hmmQueue(j) = clearHMMQueue;
                end
        end
    end

%     segmentLabel = cell(1, length(segmentStruct.segmentId));
%     for i = 1:length(segmentStruct.segmentId)
%         segmentLabel{i} = templateModel.motionName{segmentStruct.segmentId(i)};
%     end
    
    observation.currDataChecked = length(obsData)-W;

    fhmmIntermediate.stateStruct = stateStruct;
    fhmmIntermediate.clearHMMQueue = clearHMMQueue;
    fhmmIntermediate.crossingStruct = crossingStruct;
    fhmmIntermediate.veloStruct = veloStruct;
    fhmmIntermediate.segmentStruct = segmentStruct;
    
    algSegmentResult = segmentStruct;
end

function velo = veloCalc(time, data, filter, filterParam, dof)
    % expecting data to be vertical (dofs in col)
    if ~exist('filter', 'var')
        filter = 0;
    end
    
    if ~exist('dof', 'var')
        dof = size(data, 2);
    end
    
    velo = zeros(size(data));
    
    for i = 1:dof
        dataTemp = data(:, i);
        
        if ~strcmpi(filter, 'noFilter') && length(dataTemp) > 15
            dataToDiff = filterData(dataTemp, filter, filterParam); % and filter TODO
        else
            dataToDiff = dataTemp;
        end
        
        velo(:, i) = [0; diff(dataToDiff) ./ diff(time)];
    end
end

function  windowedDataSpline = resizeObservationData(windowedData, targetLength)    

    % resample the data so it's the same length as template
    dof = size(windowedData, 1);
    windowedDataLen = size(windowedData, 2);
    origSplineArray = 1:size(windowedData, 2);
    
    incriment = windowedDataLen/targetLength;
    resampleSplineArray = 1:incriment:windowedDataLen+incriment;
    windowedDataSpline = zeros(dof, length(resampleSplineArray));

    for ndof = 1:dof
        windowedDataSpline(ndof, :) = spline(origSplineArray, windowedData(ndof, :), resampleSplineArray);
    end
end

function [loglik, loglikInd, templateId] = maxLL(loglikMtxSeg, loglikMtxDof, LLThreshold, LLThresholdOverall, windowSize)

    loglikMtxNormalize = NaN(size(loglikMtxSeg));
    
    % select the max template
    [loglikMaxTemplate, loglikMaxInd] = max(loglikMtxSeg);
    templateId = loglikMtxDof(loglikMaxInd);
    
    for i = 1:length(loglikMtxSeg)
        LLVal = loglikMtxSeg(i);
        
        if loglikMtxDof(i) ~= templateId
            continue
        end
        
        % take the higher threshold
        LLThresholdTemplate = LLThreshold(loglikMtxDof(i));
        if LLThresholdTemplate > LLThresholdOverall
            threshold = LLThresholdTemplate;
        else
            threshold = LLThresholdOverall;
        end
        
        if LLVal > threshold
            % for those that are above threshold, penalize longer dataset
%             loglikMtxNormalize(i) = (LLVal - threshold) ./ windowSize(i);
            
            % don't penalize instead
            loglikMtxNormalize(i) = LLVal;
        end
    end

    [loglikNorm, loglikInd] = max(loglikMtxNormalize);

    if isempty(loglikInd) 
        % they're both empty, implying that all values were below
        % threshold, so the only thing in loglikMtxNormalize is NaNs
        [loglikNorm, loglikInd] = max(loglikMtxSeg);

        if isempty(loglikInd) 
            % none of the LL were good
            loglikInd = 1; 
        end
    end

    templateId = loglikMtxDof(loglikInd);
    loglik = loglikMtxSeg(loglikInd);
end

function [zvcIndex, zvcTime] = zvcTimeReturn_multipeak(timeData, featureCrossing, crossingMtx, crossingTime, crossingIndex, firstVelo, secondVelo, type)
    % this function is designed to be used with the peak feature extraction
    % this function will pull the first two entries of featureCrossing
    % (which it will assume is some from of ZVC crossing, followed by a
    % peak direction) and return all the ZVC values that occur before the
    % said peak. if type is set to 'last', then the reverse happens. the
    % last few ZVCs are returned
    
    % returns both the time and index values
    
    switch type
        case 'first'
            featureZVC = featureCrossing(1);
%             featurePeak = featureCrossing(2);
%             minVeloTime = max(minVeloTime);
%             maxVeloTime = max(maxVeloTime);
            timePeak = firstVelo;
            
        case 'last'
            featureZVC = featureCrossing(end);
%             featurePeak = featureCrossing(end-1);
%             minVeloTime = min(minVeloTime);
%             maxVeloTime = min(maxVeloTime);
            timePeak = secondVelo;
    end
    
%     switch featurePeak
%         case -2
%             % negative peak
%             timePeak = minVeloTime;
%             
%         case 2
%             % positive peak
%             timePeak = maxVeloTime;
%             
%         otherwise
%             timePeak = 0;
%     end

    % find all the zvc points that is the correct direction
    zvcPointsDir = find(crossingMtx == featureZVC); 
    
    switch type
        case 'first'
            % find all the entries that occurs before this time point
            zvcPointsTime = find(crossingTime <= timePeak);
            
            % include the closest 3 points to the peak, regardless of
            % sign of ZVC direction
            if length(zvcPointsTime) > 2
                zvcPointsTimeAugment = zvcPointsTime(end-2:end);
            elseif length(zvcPointsTime) > 1
                zvcPointsTimeAugment = zvcPointsTime(end-1:end);
            elseif ~isempty(zvcPointsTime)
                zvcPointsTimeAugment = zvcPointsTime(end);
            end
            
        case 'last'
            zvcPointsTime = find(crossingTime >= timePeak);
            
            if length(zvcPointsTime) > 2
                zvcPointsTimeAugment = zvcPointsTime(1:3);
            elseif length(zvcPointsTime) > 1
                zvcPointsTimeAugment = zvcPointsTime(1:2);
            elseif ~isempty(zvcPointsTime)
                zvcPointsTimeAugment = zvcPointsTime(1);
            end
    end
    
%     zvcPoints3 = unique(sort([intersect(zvcPointsDir, zvcPointsTime) zvcPointsTimeAugment])); % if there were different dir of ZVCs
    zvcPoints3 = unique(sort([zvcPointsTime zvcPointsTimeAugment])); % all the same direction of ZVC
    
    % pull out all the time points corresponding to it
    zvcIndex = crossingIndex(zvcPoints3);
    zvcTime = crossingTime(zvcPoints3);
end

%% HMM Support Functions
function [loglik, modelIndex] = hmmLogLikelihood(templateStruct, testData, currentSigDof, allActiveTemplates)
    % given a set of testData, compare against the hmModel that has been
    % previously constructed with atlasTemplateBuild and compute the
    % log-likelihood 
    
    % input parameters
    % prior1 - initial state guesses
    % transmat1 - state transition matrix, probability of
    % mu1 - mean matrix
    % sigma1 - coverance matrix (error tolerance)
    % mixmat1 - mixture matrix
    % testData - the data that will be compared against the hmModel
    
    % output
    % loglik - the log-likelihood value of the comparison
       
    % compute log-likelihood
    
%     filteredData = filterData(testData', filterFreq, filterSample, filterOrder);
%     filteredData = filteredData';
    hmmModelCount = length(templateStruct.hmmModel);
    dataLength = size(testData, 2);
    
    if ~exist('currentSigDof', 'var')
        % run forward algorithm on all HMMs
        modelIndex = 1:hmmModelCount;
    elseif exist('allActiveTemplates', 'var')
        % so allActive is passed in
        modelIndex = allActiveTemplates;
    else
        modelIndex = [];
        currentSigDof = currentSigDof(currentSigDof > 0);
        
        for i = 1:hmmModelCount
            sigdofintersect = intersect(templateStruct.featureModel{i}.sigDof, currentSigDof);
            
            if ~isempty(sigdofintersect)
                modelIndex = [modelIndex i];
            end
        end
    end
    
    loglik = zeros(1, 1, length(modelIndex));
    
    for i = 1:length(modelIndex)
%         prior =    templateStruct.hmmModel{modelIndex(i)}.prior;
%         transmat = templateStruct.hmmModel{modelIndex(i)}.transmat;
%         mu =       templateStruct.hmmModel{modelIndex(i)}.mu;
%         sigma =    templateStruct.hmmModel{modelIndex(i)}.sigma;
%         mixmat =   templateStruct.hmmModel{modelIndex(i)}.mixmat;
%         loglik(1, 1, i) = mhmm_logprob(testData, prior, transmat, mu, sigma, mixmat); 

        try
            loglik(1, 1, i) = hmmLogLik(templateStruct.hmmModel{modelIndex(i)}, testData); % or could length normalize
        catch err
            fprintf(['LL determination error, in template ' num2str(modelIndex(i)) ': ' err.message ' \n']);
            loglik(1, 1, i) = -Inf;
        end
    end
    
    modelIndex = reshape(modelIndex, 1, 1, length(modelIndex));
end

function loglik = hmmLogLik(hmmModel, testData)
    % opens the HMM structure and runs a log-likelihood comparison
    prior =    hmmModel.prior;
    transmat = hmmModel.transmat;
    mu =       hmmModel.mu;
    sigma =    hmmModel.sigma;
    mixmat =   hmmModel.mixmat;
    
    loglik = mhmm_logprob(testData, prior, transmat, mu, sigma, mixmat);
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

function veloStruct = veloCheck(timeData, veloM, t, veloStruct, crossingStruct, featureModel)
    % populate the veloStruct with peak information between ZVCs
    % zvc(1) - velo(1) - zvc(2) ... so velo will always come after a zvc,
    % and zvc array will always be one longer
    
    for i = 1:length(veloStruct)
        while veloStruct{i}.IndexZVCLastChecked < crossingStruct{i}.Counter
            
%             if veloStruct{i}.IndexZVCLastChecked == 0
%                 % hasn't been initalized yet
%                 zvcPt1 = crossingStruct{i}.Index(1);
%                 zvcPt2 = crossingStruct{i}.Index(2);
%                 
%                 veloStruct{i}.IndexZVCLastChecked = 1;
%             else
                lastCheckedInd = veloStruct{i}.IndexZVCLastChecked;
                zvcPt1 = crossingStruct{i}.Index(lastCheckedInd);
                zvcPt2 = crossingStruct{i}.Index(lastCheckedInd+1);
%             end
            
            [crossingMtx, timeLocal, veloLocal] = peakSeeker(timeData, veloM(i, :), [zvcPt1 zvcPt2]);
            
            veloStruct{i}.Index = [veloStruct{i}.Index find(timeData == timeLocal)];
            veloStruct{i}.Time = [veloStruct{i}.Time timeLocal];
            veloStruct{i}.Mtx = [veloStruct{i}.Mtx crossingMtx];
            veloStruct{i}.Mag = [veloStruct{i}.Mag veloLocal];
            veloStruct{i}.lastAttenuated = [veloStruct{i}.lastAttenuated timeLocal];
            
            veloStruct{i}.IndexZVCLastChecked = veloStruct{i}.IndexZVCLastChecked + 1;
            veloStruct{i}.Counter = length(veloStruct{i}.Time);
        end
    end
end

function veloStruct = attenuateVelo(veloStruct, currentTime, sysparam)
    % attenuates the velocity according to the sysparam
    
    for i = 1:length(veloStruct)
        for j = 1:length(veloStruct{i}.Time)
            if currentTime - veloStruct{i}.Time(j) > sysparam.feature.attenuateTime ...
                    && currentTime - veloStruct{i}.lastAttenuated(j) > sysparam.feature.attenuateRate
                
                veloStruct{i}.Mag(j) = veloStruct{i}.Mag(j) * sysparam.feature.attenuatePercent;
                veloStruct{i}.lastAttenuated(j) = currentTime;
            end
        end
    end
end

function [passFlag, peakInd, veloStruct, pivotPeakTime, pivotPeakMag] = checkVeloStruct(featureModel, crossingStruct, veloStruct, sysparam)
    % recurse through the velostructmtx TODO TODOTODOTODOTODO
    passFlag = 0;
    peakInd = [];
    pivotPeakTime = [];
    pivotPeakMag = [];
    
    % purge small magnitudes, to keep the structs small
    if featureModel.doubleZvcFlag == 0
        smallMagRejectThreshold = min(abs(featureModel.peakMag{1})) * sysparam.feature.smallMagRejectMultiplier;
    else
        smallMagRejectThreshold = min(abs(featureModel.peakMag{1}(2:end))) * sysparam.feature.smallMagRejectMultiplier;
    end
    
    veloPruneIndex = [];
    for i = 1:length(veloStruct.Mag)
        if abs(veloStruct.Mag(i)) < smallMagRejectThreshold
            veloPruneIndex = [veloPruneIndex i];
%             fprintf('Pruning small velocities \n');
            % if there's no doublezvc flag, purge the ZVCs too
            if featureModel.doubleZvcFlag == 0
                % todo figure out how to do this
            end
        end
    end
    
    veloKeepIndex = setxor(1:veloStruct.Counter, veloPruneIndex); % note the values to keep
    veloStruct.Index = veloStruct.Index(veloKeepIndex);
    veloStruct.Time = veloStruct.Time(veloKeepIndex);
    veloStruct.Mtx = veloStruct.Mtx(veloKeepIndex);
    veloStruct.Mag = veloStruct.Mag(veloKeepIndex);
    veloStruct.lastAttenuated = veloStruct.lastAttenuated(veloKeepIndex);
    veloStruct.Counter = length(veloKeepIndex);
    
    if veloStruct.Counter > featureModel.peakCount - 1
        % what if we start from the start of an unchecked cycle?
        veloStruct.IndexVeloLastChecked = 1;
        
        currInd = veloStruct.IndexVeloLastChecked;
        currPeak = 0;
        
        [veloStruct, peakInd, passFlag, pivotPeakTime, pivotPeakMag] = checkVeloRecurse(veloStruct, crossingStruct, featureModel, currInd, currPeak, sysparam);
    end
end

function [veloStruct, peakInd, passFlag, pivotPeakTime, pivotPeakMag] = checkVeloRecurse(veloStruct, crossingStruct, featureModel, currInd, currPeak, sysparam)

    featurePeakCount = featureModel.peakCount;
    peakFeature = featureModel.crossingMtx{1};
    peakTimeThresholdMultiplier = sysparam.feature.peaktopeakMultiplier;
    
%     smallMagMultiplier = sysparam.feature.velocityMultiplierThreshold; 
    peakToPeakTimeThreshold = (featureModel.veloTime{1}(end)  - featureModel.veloTime{1}(1)) * sysparam.feature.peaktopeakMultiplier; %TODO

    pivotPeakTime = [];
    pivotPeakMag = [];
    
    if ~exist('peakInd', 'var')
        peakInd = [];
        cumCurrPeak = [];
    end
    
    passFlag = 0;
    lastInd1 = currInd; % initalizing values
    lastInd2 = lastInd1 + 1;
    peakInd = zeros(1, featurePeakCount);
    
    ind = find(abs(veloStruct.Mag) > 0);
    veloStruct.Index = veloStruct.Index(ind);
    veloStruct.Time = veloStruct.Time(ind);
    veloStruct.Mtx = veloStruct.Mtx(ind);
    veloStruct.Mag = veloStruct.Mag(ind);
    veloStruct.lastAttenuated = veloStruct.lastAttenuated(ind);
    veloStruct.Counter = length(ind);
    
    % combine the velocity and feature matrix
    combinedFeatureIndex = [crossingStruct.Index veloStruct.Index];
    combinedFeatureTime = [crossingStruct.Time veloStruct.Time];
    combinedFeatureMtx = [crossingStruct.Mtx veloStruct.Mtx];
    combinedFeatureMag = [zeros(size(crossingStruct.Mtx)) veloStruct.Mag];
    
    [detectedFeatureTime, detectedFeatureInd] = sort(combinedFeatureTime);
    detectedFeatureIndex = combinedFeatureIndex(detectedFeatureInd);
    detectedFeatureMtx = combinedFeatureMtx(detectedFeatureInd);
    detectedFeatureMag = combinedFeatureMag(detectedFeatureInd);
    detectedFeatureCounter = length(detectedFeatureTime);
    
    passFlag = 0;

%     featureItemToMatch = featureModel.crossingMtx{1};
%     
%     % actually, lets try this
%     for i = length(featureItemToMatch):length(detectedFeatureMtx)
%         currentFeatureSet = detectedFeatureMtx(i-length(featureItemToMatch)+1:i);
%         
%         if isequal(currentFeatureSet, featureItemToMatch)
%             % found a matching sequence
%             peakInd = [];
%             
%             matchedVelo = detectedFeatureMag(i-length(featureItemToMatch)+1:i);
%             nonZeroVelo = matchedVelo(matchedVelo ~= 0);
%             
%             pivotPeak = [nonZeroVelo(1) nonZeroVelo(end)];
%             passFlag = 1;
%             return
%         end
%     end
    
    
    crossingMtxInd = 0; % 
    crossingMtxMatch = []; % record the time index of where the matches occured
    lastIndexExamined = 0; % position in the detectedFeatureIndex array
    crossingMtxInd = crossingMtxInd + 1;
    
    % okay. lets try this
    % search for the peaks
    featureModelVeloInd = find(featureModel.crossingMtx{1} ~= 1);
    featurePeakItems = featureModel.crossingMtx{1}(featureModelVeloInd); % pull out all non-ZVC entries
    veloMtxInd = 0;
    passVeloCounter = 0;
    passVeloIndArray = {}; % the ind of the pass velos
    
    for i = length(featurePeakItems):length(veloStruct.Mtx)
        currIndex = i-length(featurePeakItems)+1:i;
        currentFeatureSet = veloStruct.Mtx(currIndex);
        
        if isequal(currentFeatureSet, featurePeakItems)
            % found a matching sequence
            peakVeloInd = veloStruct.Index(currIndex);
            peakVeloTime = veloStruct.Time(currIndex);
            peakVeloMag = veloStruct.Mag(currIndex);

            % search the front
            currentFeatureSet = featureModel.crossingMtx{1}(1:featureModelVeloInd(1)-1);
            currentSubset = crossingStruct.Index(crossingStruct.Index <= peakVeloInd(1));
            if length(currentSubset) >= length(currentFeatureSet)
                startIndTmp = find(currentSubset, length(currentFeatureSet), 'last');
                startInd = currentSubset(startIndTmp(1));
            else
                % don't have enough ZVCs...
                fprintf('Not enough ZVCs with properly sized velocities\n');
                continue
            end
            
            % search the end
            currentFeatureSet = featureModel.crossingMtx{1}(featureModelVeloInd(end)+1:end);
            currentSubset = crossingStruct.Index(crossingStruct.Index >= peakVeloInd(end));
            if length(currentSubset) >= length(currentFeatureSet)
                endIndTmp = find(currentSubset, length(currentFeatureSet), 'first');
                endInd = currentSubset(endIndTmp(1));
            else
                fprintf('Not enough ZVCs with properly sized velocities\n');
                continue
            end
            
            % make the sure the peak to peak time (which should be first
            % peak to last peak) is appropriate
            if peakVeloTime(end) - peakVeloTime(1) < peakToPeakTimeThreshold
                % smaller then the threshold
                fprintf('Peak-to-peak time too small\n');
                continue
            end

            passVeloCounter = passVeloCounter + 1;
            passVeloIndArray{passVeloCounter, 1} = peakVeloTime; % time to the velocities
            passVeloIndArray{passVeloCounter, 2} = peakVeloMag; % total abs largest, used to compare between them           
        end
    end
    
%     if passVeloCounter % found velocity matches 
%         [abssumVeloVal, absSumVeloInd] = max(cell2mat(passVeloIndArray(:, 2))); % pull the fastest one
%         peakInd = [];
%         pivotPeak = passVeloIndArray(absSumVeloInd, 1);
%         passFlag = 1;        
%     end
    
    if passVeloCounter % found velocity matches
        peakInd = [];
        pivotPeakTime = passVeloIndArray(:, 1);
        pivotPeakMag = passVeloIndArray(:, 2);
        passFlag = 1;        
    end

end

function [crossingStruct, veloStruct] = cleanCrossingStructs(crossingStruct, veloStruct, recycleTime)
    indexToKeep = find(crossingStruct.Time >= recycleTime);
    crossingStruct.Index = crossingStruct.Index(indexToKeep);
    crossingStruct.Time = crossingStruct.Time(indexToKeep);
    crossingStruct.Mtx = crossingStruct.Mtx(indexToKeep);
    crossingStruct.Counter = length(indexToKeep);

    if ~isempty(crossingStruct.Time)
        veloKeepTime = crossingStruct.Time(1);

        indexToKeep = find(veloStruct.Time >= veloKeepTime);
        veloStruct.Index = veloStruct.Index(indexToKeep);
        veloStruct.Time = veloStruct.Time(indexToKeep);
        veloStruct.Mtx = veloStruct.Mtx(indexToKeep);
        veloStruct.Mag = veloStruct.Mag(indexToKeep);
        veloStruct.lastAttenuated = veloStruct.lastAttenuated(indexToKeep);
        veloStruct.Counter = length(indexToKeep);
    else
        veloStruct.Index = [];
        veloStruct.Time = [];
        veloStruct.Mtx = [];
        veloStruct.Mag = [];
        veloStruct.lastAttenuated = [];
        veloStruct.Counter = length(indexToKeep);
    end


    veloStruct.IndexVeloPass = [];
    if crossingStruct.Counter > 0
        veloStruct.IndexZVCLastChecked = crossingStruct.Counter;
    else
        veloStruct.IndexZVCLastChecked = 1;
    end

    if veloStruct.Counter > 0
        veloStruct.IndexVeloLastChecked = veloStruct.Counter;
    else
        veloStruct.IndexVeloLastChecked = 1;
    end
end