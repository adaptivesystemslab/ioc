function [jointAngle, jointTime, segmentInfo] = dataSource_TRO2009(dataSettings, pathToRawData, motionSet, Ntemplate)
    % load and process Dana's TRO dataset (long)
    % remodified to fit new featurem-hmm format
    
    if ~exist('hmm', 'var')
        hmm = 1;
    end
    
    Nsegmentpass = 345; % for observation data, pass this number segments (of 349)
%     Ntemplate = 10; % number of templates of set (usable template)
    frameRate = 33 + 1/3;
%     frameRate = 60;
    combineFlexExtMotions = 1;
    
    dofToUse = 9:28; % full
%     dofToUse = [21:28]; % right arm and left arm 21 22 25 26
%     dofToUse = [9 18 21 22 25 26]; % right arm and left arm and upper leg markers
    
    %load joint data
    baseFolder = pathToRawData; 
    load(fullfile(baseFolder, 'tro2009_joint.mat'));
    
    % load segment data
%     segmentData = motionAnalysis_TRO2009_segFileRead;
    [segmentBreakdown, segmentSorted] = segmentAnalysis(segmentData.segment(:, 3));
    
    % copy out joint names
    jointNames = jointData.marker(dofToUse);
    
    templateData = cell(1, length(motionSet));
%     templateTime = cell(1, length(motionSet));

    segmentInd = cell2mat(segmentData.segment(1:end-2, 1:2));
    segmentLabels = segmentData.segment(1:end-2, 3);
    segmentId = cell2mat(segmentData.segment(1:end-2, 4));
    
    segmentStart = [];
    segmentEnd = [];
    segmentName = {};
    
    segmentStartStartToEnd = [];
    segmentEndStartToEnd = [];
    
    % combine the templates
    for i = 1:length(motionSet)
        % pull, filter
%         templateData{i} = dataSetTemplate(jointData.joint, segmentInd, segmentLabels, segmentId, motionSet{i}, i, Ntemplate, frameRate, dofToUse, hmm);
        [timeStart, timeEnd, segmentNameTemp, timeStartStartToEnd, timeEndStartToEnd] = dataSetTemplate_pointsOnly(jointData.joint(:, dofToUse), segmentLabels, segmentInd, segmentId, motionSet{i}, Ntemplate, combineFlexExtMotions);
        segmentStart = [segmentStart; timeStart];
        segmentEnd   = [segmentEnd;   timeEnd];
        segmentName  = [segmentName; segmentNameTemp];
        
        segmentStartStartToEnd = [segmentStartStartToEnd; timeStartStartToEnd];
        segmentEndStartToEnd   = [segmentEndStartToEnd;   timeEndStartToEnd];
    end
    
    [c, ind] = sort(segmentStart);
    
    segmentInfo.timeStart = segmentStart(ind)/frameRate;
    segmentInfo.timeEnd = segmentEnd(ind)/frameRate;
    segmentInfo.segmentName = segmentName(ind);
    segmentInfo.timeStartOneBefore = segmentStartStartToEnd;
    segmentInfo.timeEndOneAfter = segmentEndStartToEnd;
    
%     debug = [segmentInfo.TimeStart segmentInfo.TimeEnd segmentInfo.TimeEnd - segmentInfo.TimeStart];
    
    % pull out a certain amount of segments 
%     segIndSegStart = segmentStart(1) - 50; % start from the first segment
%     segIndSegEnd = segmentEnd(end) + 50;

    segIndSegStart = 1; % start from the first segment
    segIndSegEnd = length(jointData.joint);
    
    jointTime = (segIndSegStart:segIndSegEnd)'/frameRate;
    jointAngle = unwrap(jointData.joint(segIndSegStart:segIndSegEnd, dofToUse)); % stripping frame count
    
    
    
%     obsTime{1} = obsTimeRaw;
%     obsData{1} = obsDataRaw;
%     
%     segmentIdArray = [];
%     segmentTimeArray = [];
%     segmentLabelArray = {};
%     segmentCountTotal = 0;
%     
%     % pull the manual segmentation set for this
%     for i = 1:length(motionSet)
%         % pull, filter
%         [segmentTime, segmentLabel, segmentCount] = dataSetSegmentation(segmentInd, segmentLabels, segmentId, motionSet{i}, Nsegmentpass, frameRate);
%         
%         segmentTimeArray = [segmentTimeArray segmentTime];
%         segmentCountTotal = segmentCountTotal + segmentCount;
%         segmentLabelArray = [segmentLabelArray segmentLabel];
%         segmentIdArray = [segmentIdArray i*ones(size(segmentLabel))];
%     end
%     
%     [segmentTimeArray, segmentTimeMapping] = sort(segmentTimeArray);
%     
% %     manualSeg{1}.Time = segmentTimeArray;
% %     manualSeg{1}.Count = segmentCountTotal;
% %     manualSeg{1}.Label = segmentLabelArray(segmentTimeMapping);
% %     manualSeg{1}.Id = segmentIdArray(segmentTimeMapping);
%     
%     observationData{1}.data = obsData;
%     observationData{1}.time = obsTime;
%     observationData{1}.label = segmentLabelArray(segmentTimeMapping);
%     observationData{1}.id = segmentIdArray(segmentTimeMapping);
%     observationData{1}.manualSegTime = segmentTimeArray;
%     observationData{1}.count = segmentCountTotal;
end

function [timeStart, timeEnd, segmentName, timeStartStartToEnd, timeEndStartToEnd] = dataSetTemplate_pointsOnly(dataset, segmentLabels, segmentInd, segmentId, motionName, templateCount, combineFlexExtMotions)
    % pull out the template motion, the last 'templateCount' few

    % BW filter parameters
%     filterFreq = 0.5;
%     filterSample = 25;
%     filterOrder = 5;
    
    [motionLocations, motionLocationsStartToEnd, motionId, motionNameArray] = motionLocator(segmentLabels, segmentId, motionName, combineFlexExtMotions);
    
    if size(motionLocations, 1) < templateCount
        fprintf(['DataSource2: Template available less than desired for ' ...
            motionName, ', have ' num2str(size(motionLocations, 1)) ' but want ' num2str(templateCount) '\n']);
        templateCount = size(motionLocations, 1);
    end
    
    motionLocationCount = motionLocations(end-templateCount+1:end, :);
    timeStart = zeros(size(motionLocationCount, 1), 1);
    timeEnd = zeros(size(motionLocationCount, 1), 1);
    segmentName = motionNameArray(end-templateCount+1:end, :);
    
    for i = 1:size(motionLocationCount, 1)
        timeStart(i) = segmentInd(motionLocationCount(i, 1), 1);
        timeEnd(i)   = segmentInd(motionLocationCount(i, end), 2);
%         segmentName{i} = motionName{i};
    end
    
    motionLocationsStartToEndCount = motionLocationsStartToEnd(end-templateCount+1:end, :);
    timeStartStartToEnd = zeros(size(motionLocationCount, 1), 1);
    timeEndStartToEnd = zeros(size(motionLocationCount, 1), 1);
    
    for i = 1:size(motionLocationCount, 1)
        timeStartStartToEnd(i) = segmentInd(motionLocationsStartToEndCount(i, 1), 1);
        timeEndStartToEnd(i)   = segmentInd(motionLocationsStartToEndCount(i, end), 2);
%         segmentName{i} = motionName{i};
    end
    
    templateData = unwrap(dataset(timeStart(1):timeEnd(1), :)); % strip out the frame number, and unwrap angle
    
%     templateTime = cell(1, templateCount);
%     templateDataSet = cell(1, templateCount);
%     segTime = zeros(1, templateCount*2);
%     label = cell(1, templateCount*2);
%     id = zeros(1, templateCount*2);
%     
%     for i = 1:templateCount
%         segStart = segmentInd(motionLocationCount(i, 1), 1); % combine the "up" and "down" motion
%         segEnd = segmentInd(motionLocationCount(i, end), 2);
%         time = (segStart:segEnd)'/frameRate;
%         
%         templateData = unwrap(dataset(segStart:segEnd, dofToUse)); % strip out the frame number, and unwrap angle        
% %         templateData = filter_dualpassBW(templateData, filterFreq, filterSample, filterOrder); % and filter
%         
%         if hmm
%             time = time';
%             templateData = templateData';
%         end
%         
%         templateTime{i} = time;
%         templateDataSet{i} = templateData;
%         
%         doubleInd = (i-1)*2+1:(i)*2;
%         segTime(doubleInd) = time([1 end]);
%         label{doubleInd(1)} = motionName;
%         label{doubleInd(2)} = motionName;
%         id(doubleInd) = [motionNameId motionNameId];
% 
%     end
%     
%     templateStruct.data = templateDataSet;
%     templateStruct.time = templateTime;
%     templateStruct.label = label;
%     templateStruct.id = id;
%     templateStruct.manualSegTime = segTime;
%     templateStruct.count = templateCount;
end
% 
% function templateStruct = dataSetTemplate(dataset, segmentInd, segmentLabels, segmentId, motionName, motionNameId, templateCount, frameRate, dofToUse, hmm)
%     % pull out the template motion, the last 'templateCount' few
% 
%     % BW filter parameters
% %     filterFreq = 0.5;
% %     filterSample = 25;
% %     filterOrder = 5;
%     
%     [motionLocations, motionId, segmentName] = motionLocator(segmentLabels, segmentId, motionName);
%     
%     if size(motionLocations, 1) < templateCount
%         fprintf(['DataSource2: Template available less than desired for ' ...
%             motionName, ', have ' num2str(size(motionLocations, 1)) ' but want ' num2str(templateCount) '\n']);
%         templateCount = size(motionLocations, 1);
%     end
%     
%     motionLocationCount = motionLocations(end-templateCount+1:end, :);
%     
%     templateTime = cell(1, templateCount);
%     templateDataSet = cell(1, templateCount);
%     segTime = zeros(1, templateCount*2);
%     label = cell(1, templateCount*2);
%     id = zeros(1, templateCount*2);
%     
%     for i = 1:templateCount
%         segStart = segmentInd(motionLocationCount(i, 1), 1); % combine the "up" and "down" motion
%         segEnd = segmentInd(motionLocationCount(i, end), 2);
%         time = (segStart:segEnd)'/frameRate;
%         
%         templateData = unwrap(dataset(segStart:segEnd, dofToUse)); % strip out the frame number, and unwrap angle        
% %         templateData = filter_dualpassBW(templateData, filterFreq, filterSample, filterOrder); % and filter
%         
%         if hmm
%             time = time';
%             templateData = templateData';
%         end
%         
%         templateTime{i} = time;
%         templateDataSet{i} = templateData;
%         
%         doubleInd = (i-1)*2+1:(i)*2;
%         segTime(doubleInd) = time([1 end]);
%         label{doubleInd(1)} = motionName;
%         label{doubleInd(2)} = motionName;
%         id(doubleInd) = [motionNameId motionNameId];
% 
%     end
%     
%     templateStruct.data = templateDataSet;
%     templateStruct.time = templateTime;
%     templateStruct.label = label;
%     templateStruct.id = id;
%     templateStruct.manualSegTime = segTime;
%     templateStruct.count = templateCount;
% end
% 
% function [segmentTime, segmentLabel, segmentCount] = dataSetSegmentation(segmentInd, segmentLabels, segmentId, motionName, segmentPassCount, frameRate)
%     % pull out the template motion, the last 'templateCount' few
% 
% %     % BW filter parameters
% %     filterFreq = 0.5;
% %     filterSample = 25;
% %     filterOrder = 5;
%     
%     [motionLocations, namedSegmentId, segmentName] = motionLocator(segmentLabels, segmentId, motionName);
%     motionLocationInd = find(namedSegmentId(:, 1) <= segmentPassCount);
%     motionLocationCount = motionLocations(motionLocationInd, :);
%     segmentCount = length(motionLocationInd);
%     
%     segmentTime = zeros(1, segmentCount*2);
%     segmentLabel = cell(1, segmentCount*2);
%     
%     for i = 1:length(motionLocationInd)
%         segStart = segmentInd(motionLocationCount(i, 1), 1); % combine the "up" and "down" motion
%         segEnd = segmentInd(motionLocationCount(i, end), 2);
% %         time = (segStart:segEnd)'/frameRate;
%         
%         segmentTime(2*i - 1) = segStart/frameRate;
%         segmentTime(2*i) = segEnd/frameRate;
%         
%         segmentLabel{2*i - 1} = motionName;
%         segmentLabel{2*i} = motionName;
%     end
% end
% 
function [motionLocations, motionLocationsStartToEnd, segmentId, segmentName] = motionLocator(segmentLabels, segmentCount, desiredSegment, combineFlexExtMotions)
    % find all the occurances of the 'desiredSegment'
    motionLocations = [];
    segmentId = [];
    segmentCountExpended = []; % once the segment has been consumed, add it here
    motionLocationsStartToEnd = [];
    segmentIndStartToEnd = [];
    
    for i = 1:length(segmentLabels)
% % %         if strcmpi(segmentLabels{i}, desiredSegment)
% % % %             segmentInd = find(segmentCount == segmentCount(i)); % pull all entries that match that ID
% % % %             endLoc = find(segmentCount == segmentCount(i), 2, 'last'); % pull just the last 2
% % %             segmentInd = [i; i+1];
% % %             
% % %             if isempty(motionLocations)
% % %                 motionLocations = segmentInd';
% % %                 segmentId = segmentCount(i)*ones(size(segmentInd))';
% % %             else
% % % %                 if size(motionLocations, 2) ~= size(segmentInd', 2)
% % % %                     segmentInd = [segmentInd; 0];
% % % %                 end
% % %                 
% % % 
% % %                 motionLocations = [motionLocations; segmentInd'];
% % %                 segmentId = [segmentId; segmentCount(i)*ones(size(segmentInd))'];
% % %             end

        if strcmpi(segmentLabels{i}, desiredSegment) || ...
                (strcmpi('ALLMOTIONS', desiredSegment) && ~sum(segmentCountExpended == segmentCount(i)))
            currSegmentCount = segmentCount(i);
            if combineFlexExtMotions
                segCountInd = find(segmentCount == currSegmentCount);
                segmentInd = [segCountInd(1); segCountInd(end)];
                segmentIndStartToEnd = [segCountInd(1)-1; segCountInd(end)+1];
            else
                segmentInd = [i; i];
            end
            segmentCountExpended = [segmentCountExpended, currSegmentCount];
            
            if isempty(motionLocations)
                motionLocations = segmentInd';
                motionLocationsStartToEnd = segmentIndStartToEnd';
                segmentId = segmentCount(i)*ones(size(segmentInd))';
                segmentName{1} = segmentLabels{i};
            else
                motionLocations = [motionLocations; segmentInd'];
                motionLocationsStartToEnd = [motionLocationsStartToEnd; segmentIndStartToEnd'];
                segmentId = [segmentId; segmentCount(i)*ones(size(segmentInd))'];
                segmentName{end+1} = segmentLabels{i};
            end
        end
    end
    
    segmentName = segmentName';
end

function [segmentBreakdown, segmentSorted] = segmentAnalysis(segmentNames)
    segSorted = sort(segmentNames);
    segUnique = unique(segSorted, 'first');
    segCount = zeros(size(segUnique));
    
    for i = 1:length(segSorted)
        ind = find(strcmpi(segSorted{i}, segUnique));
        segCount(ind) = segCount(ind) + 1;
    end
    
    [sortedVal, sortMap] = sort(segCount, 1, 'descend');
    sortSegUnique = segUnique(sortMap);
    
    segmentBreakdown = cell(length(segUnique), 2);
    segmentSorted = cell(length(segUnique), 2);
    
    segmentBreakdown(:, 1) = segUnique;
    segmentSorted(:, 1) = sortSegUnique;
    
    for i = 1:length(segUnique)
        segmentBreakdown{i, 2} = segCount(i);
        segmentSorted{i, 2} = sortedVal(i);
    end
end