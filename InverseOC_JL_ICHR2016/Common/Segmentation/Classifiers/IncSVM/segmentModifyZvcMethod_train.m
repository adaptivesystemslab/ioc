function [segmentInfo, templateModel, runningTemplateData] = segmentModifyZvcMethod_train(mansegLoadMode, observationNew, inputSegmentInfo, model, runningTemplateData)
% now, if the dataset is training, we can replace it with an
% semi-supervised method, but the 'testing' should remain
% consistent.

if ~isempty(observationNew)
    time = observationNew.obsTime;
    [exemplar, fhmmIntermediate, sysparam] = fgInitialize;
    [exemplar, fhmmIntermediate] = updateDataStruct(exemplar, observationNew, sysparam, fhmmIntermediate);
end

if ~isempty(inputSegmentInfo)
    [~, startInd] = findClosestValue(inputSegmentInfo.timeStart, time);
    [~, endInd] = findClosestValue(inputSegmentInfo.timeEnd, time);
end


%                             for indVelo = 1:size(templateVelocity, 1)
%                                 templateVeloLocal(indVelo) = norm(templateVelocity(indVelo, sigDof));
%                             end

switch mansegLoadMode
    case 'loadAll'
        % normal. load all data available
        segmentInfo = inputSegmentInfo;
        
        trainingInd = 1:length(startInd);
        for i = trainingInd
            trainingTime{i} = exemplar.obsTime(startInd(i):endInd(i))' - time(startInd(i));
            trainingVelo{i} = exemplar.obsVelo(:, startInd(i):endInd(i));
            trainingLabel{i} = inputSegmentInfo.segmentName{i};
        end
        
        % copy them out to build fgModel 
        runningTemplateData.trainingTime = [runningTemplateData.trainingTime trainingTime];
        runningTemplateData.trainingVelo = [runningTemplateData.trainingVelo trainingVelo];
        runningTemplateData.trainingLabel = [runningTemplateData.trainingLabel trainingLabel];

        templateModel = [];
        
    case {'loadFirst1', 'loadFirst3', 'loadFirst5'}
        % semi-supervised. load the first segment, then try to
        % learn from the data. use the zvc fg model
        
        counter = 0;
        % how much data to load into training
%         trainingInd = 1:3;
%         for i = trainingInd
%             counter = counter + 1;
%             trainingTime{i} = exemplar.obsTime(startInd(i):endInd(i))' - time(startInd(i));
%             trainingVelo{i} = exemplar.obsVelo(:, startInd(i):endInd(i));
%             trainingLabel{i} = 'KEFO_SIT';
%         end
        
        switch mansegLoadMode
            case 'loadFirst1'
                trainingInd1 = [1];
                trainingInd2 = [2];
                
            case 'loadFirst3'
                trainingInd1 = [1 3 5];
                trainingInd2 = [2 4 6];
                
            case 'loadFirst5'
                trainingInd1 = [1 3 5 7 9];
                trainingInd2 = [2 4 6 8 10];
        end

        for i = trainingInd1
            counter = counter + 1;
            trainingTime{i} = exemplar.obsTime(startInd(i):endInd(i))' - time(startInd(i));
            trainingVelo{i} = exemplar.obsVelo(:, startInd(i):endInd(i));
%             trainingLabel{i} = [observationNew.obsName '_STARTING'];
            trainingLabel{i} = inputSegmentInfo.segmentName{i};
        end
        
        for i = trainingInd2
            counter = counter + 1;
            trainingTime{i} = exemplar.obsTime(startInd(i):endInd(i))' - time(startInd(i));
            trainingVelo{i} = exemplar.obsVelo(:, startInd(i):endInd(i));
%             trainingLabel{i} = [observationNew.obsName '_ENDING'];
            trainingLabel{i} = inputSegmentInfo.segmentName{i};
        end
        
        trainingInd = sort([trainingInd1 trainingInd2]);
        
        runningTemplateData.trainingTime = [runningTemplateData.trainingTime trainingTime];
        runningTemplateData.trainingVelo = [runningTemplateData.trainingVelo trainingVelo];
        runningTemplateData.trainingLabel = [runningTemplateData.trainingLabel trainingLabel];

%         featureModel = crossingPointDetection_combined('', trainingTime, trainingVelo, 1, 1, 'twoPeak');
        featureModel = crossingPointDetection_combined(unique(trainingLabel), trainingTime, trainingVelo, trainingLabel, 'multiPeak');
        templateModel = setupForFg(featureModel);
        
        % now pull out all the segments
%         obsStartInd = floor(mean([endInd(1) startInd(2)]));
        obsStartDiff = diff([startInd(trainingInd(end)) endInd(trainingInd(end))]);
        obsStartInd = floor(mean([startInd(trainingInd(end)) endInd(trainingInd(end))]) + obsStartDiff/4);
        obsEndInd = length(time);
        
        observation.obsTime = exemplar.obsTime(obsStartInd:obsEndInd) - time(obsStartInd);
        observation.obsData = exemplar.obsData(:, obsStartInd:obsEndInd);
        observation.obsVelo = exemplar.obsVelo(:, obsStartInd:obsEndInd);
        observation.currDataLength = length(observation.obsTime);
        observation.currDataChecked = 0;
        
        [segResults, ~, ~, ~] = fgSegment(observation, templateModel, sysparam, fhmmIntermediate);
        
        % add in the ground truth point to the points obtained from the fg
        % system
        startInd = [startInd(trainingInd); segResults(:, 1) + obsStartInd];
        endInd = [endInd(trainingInd); segResults(:, 2) + obsStartInd];
        
        % readd the obstime init offset back in TODO
        
%         if length(featureModel.zvcTime{1}) < 2
%             % if we don't end up with enough matches, enforce a 2 peak
%             % approach
%             
%         end
        
        % look at the crossing struct
%         switch windowMode
%             case 1
%                 % combine flex with ext
%                 [startVal, startInd, ~, startDiff] = findClosestValue(featureModel.zvcTime{1}(1)*1000, jointTime);
%                 [endVal, endInd, ~, endDiff]  = findClosestValue(featureModel.zvcTime{1}(4)*1000, jointTime);
%                 
%             otherwise
%                 % separate flex from ext
%                 [startVal, startInd(1), ~, startDiff] = findClosestValue(featureModel.zvcTime{1}(1)*1000, jointTime);
%                 [endVal, endInd(1), ~, endDiff]  = findClosestValue(featureModel.zvcTime{1}(2)*1000, jointTime);
%                 
%                 [startVal, startInd(2), ~, startDiff] = findClosestValue(featureModel.zvcTime{1}(3)*1000, jointTime);
%                 [endVal, endInd(2), ~, endDiff]  = findClosestValue(featureModel.zvcTime{1}(4)*1000, jointTime);
%         end
        
%         [startInd, endInd] = checkStartEndInd(startInd, endInd, length(time));
        
        segmentInfo.segmentCount = 1:length(startInd);
        segmentInfo.use = ones(size(startInd));
        segmentInfo.timeStart = time(startInd);
        segmentInfo.timeEnd = time(endInd);
        
        % now apply the model to look at the rest of the dataset and
        % generate segment info
        
    case 'rehashModel' % retrain all the runningTemplateData into a model
%         featureModel = crossingPointDetection_combined('', trainingTime, trainingVelo, 1, 1, 'twoPeak');
        trainingLabelUnique = unique(runningTemplateData.trainingLabel);
        featureModel = crossingPointDetection_combined(trainingLabelUnique, runningTemplateData.trainingTime, runningTemplateData.trainingVelo, runningTemplateData.trainingLabel, 'multiPeak');
        %                             [crossingStruct, crossingDof] = zvcCheck(jointTime, normVelo);
        
        templateModel = setupForFg(featureModel);
        segmentInfo = [];
end

if 0
    initTime = time(1);
    pos = exemplar.obsData';
    velo = exemplar.obsVelo';
    
    figure;
    h(1) = subplot(211);
    plot(time - initTime, pos);
    plotBoxes(h(1), inputSegmentInfo.timeStart  - initTime,inputSegmentInfo.timeEnd  - initTime, 'r', 0);
    plotBoxes(h(1), segmentInfo.timeStart  - initTime,segmentInfo.timeEnd  - initTime, 'b', 0.2);
    
    h(2) = subplot(212);
    plot(time  - initTime, velo); hold on
    plot(time  - initTime, zeros(size(time)), 'k');
    plotBoxes(h(2), inputSegmentInfo.timeStart  - initTime,inputSegmentInfo.timeEnd  - initTime, 'r', 0);
    plotBoxes(h(2), segmentInfo.timeStart  - initTime,segmentInfo.timeEnd  - initTime, 'b', 0.02);
    
    linkaxes(h, 'x');
end

% segmentInfo.SegmentCount = segmentInfo.SegmentCount(cropStart:cropEnd);
% segmentInfo.Use = segmentInfo.Use(cropStart:cropEnd);
% segmentInfo.timeStart = segmentInfo.timeStart(cropStart:cropEnd);
% segmentInfo.timeEnd = segmentInfo.timeEnd(cropStart:cropEnd);

end