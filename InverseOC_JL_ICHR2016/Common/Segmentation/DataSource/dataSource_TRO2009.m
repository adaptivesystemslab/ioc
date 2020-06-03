function [jointAngle, jointTime, segmentInfo] = dataSource_TRO2009(dataSettings, pathToRawData, motionSet, Ntemplate)
    % load and process Dana's TRO dataset (long)
    % remodified to fit new featurem-hmm format
    % segmentIncludeInd is '1' for include (train using this)
    %                      '0' for meh (don't train, but can test)
    %                     '-1' for blacklist (don't train or test)
    
    fullSegmentsFlag = dataSettings.fullSegmentFlag;
    
%     Nsegmentpass = 345; % for observation data, pass this number segments (of 349)
%     Ntemplate = 10; % number of templates of set (usable template)
    frameRate = 33 + 1/3;
%     frameRate = 60;
    
    dofToUse = 9:28; % full
%     dofToUse = [21:28]; % right arm and left arm 21 22 25 26
%     dofToUse = [9 18 21 22 25 26]; % right arm and left arm and upper leg markers

% dofToUse = [9 15 21 25]; dofToUse = [9 15 21 23 25 27];
%  9: RightLeg1
% 15: LeftLeg1
% 18: LLJ5
% 21: RightArm1
% 22: RAJ3
% 23: RightArm2
% 25: LeftArm1
% 26: LAJ3
% 27: LeftArm2

    %load joint data
    baseFolder = pathToRawData; % tro2009_joint, tro2009_joint_errorsRemoved
    jointPath = fullfile(baseFolder, 'tro2009_joint.mat'); % get segmentData and jointData, use 'motionAnalysis_TRO2009_segFileRead' to regen
    load(jointPath);
    
%     segmentInd = cell2mat(segmentData.segment(1:end, 1:2));
%     segmentLabels = segmentData.segment(1:end, 3);
%     segmentId = cell2mat(segmentData.segment(1:end, 4));
%         templateBlacklist = [];

% % %     % load segment data
    segPath = fullfile(baseFolder, 'KulicEtAlTRO_Segmentation_ManuallyLabeled_consolidated.csv')
    segmentData = motionAnalysis_TRO2009_segFileRead(segPath); % though, there is an instance attached to jointPath
%     save(fullfile(baseFolder, 'tro2009_joint.mat'), segmentData);
    
    % copy out joint names
    jointNames = jointData.marker(dofToUse);

    segmentInd = segmentData.segmentTimeModified; % segmentTimeOrig segmentData.segmentTimeModified
    segmentLabels = segmentData.segmentName;
    segmentId = segmentData.segmentId;
    
    if fullSegmentsFlag
        % full segment, blacklist by id
        templateBlacklist = [2, 3, 15, 34, 74, 181, 208, 230, 233]; % motion blacklisted due to manual seg error
    else
        % full segment, blacklist by ind
%         segmentData.segmentBlackListInd
        templateBlacklist = segmentData.segmentBlackListInd;
    end
    
    

    [segmentBreakdown, segmentSorted] = segmentAnalysis(segmentLabels);
    
%     h = figure;
%     plot(jointData.joint(:, dofToUse));
%     plotBoxes(h, segmentInd(:, 1), segmentInd(:, 2));
    
    
%     % combine the templates    
    [segmentStart, segmentEnd, segmentName, segmentIncludeInd] = dataSetTemplate_pointsOnly(jointData.joint(:, dofToUse), segmentLabels, segmentInd, segmentId, motionSet, Ntemplate, templateBlacklist, fullSegmentsFlag);

    % incorporate blacklist
%     segmentIncludeInd = +segmentIncludeInd; % remove the logic attribute
%     segmentIncludeInd(1) = segmentIncludeInd(1) - 1;
    segmentIncludeInd(templateBlacklist) = -1;
    
    segmentInfo.timeStart = segmentStart/frameRate;
    segmentInfo.timeEnd = segmentEnd/frameRate;
    segmentInfo.segmentName = segmentName;
    segmentInfo.segmentIncludeInd = segmentIncludeInd; 
    
    
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

function [timeStart, timeEnd, segmentName, includeVector] = dataSetTemplate_pointsOnly(dataset, segmentLabels, segmentInd, segmentId, motionName, templateCount, templateBlacklist, fullSegmentsFlag)
    % pull out the template motion, the last 'templateCount' few

    % BW filter parameters
%     filterFreq = 0.5;
%     filterSample = 25;
%     filterOrder = 5;
    
    % pull out all the motions and identify the timestamps
    [motionLocationsFull, motionIdFull, segmentName] = motionLocator(segmentLabels, segmentId, 'ALLMOTIONS', fullSegmentsFlag);
    
    timeStart = zeros(size(motionLocationsFull, 1), 1);
    timeEnd = zeros(size(motionLocationsFull, 1), 1);
    
    for i = 1:size(motionLocationsFull, 1)
        timeStart(i) = segmentInd(motionLocationsFull(i, 1), 1);
        timeEnd(i)   = segmentInd(motionLocationsFull(i, end), 2);
    end

    includeMatrix = zeros(length(motionLocationsFull), length(motionName));
    for i = 1:length(motionName)
        currMotion = motionName{i};
        [motionLocations, motionId, motionNameArray] = motionLocator(segmentLabels, segmentId, currMotion, fullSegmentsFlag);
        
        
        if fullSegmentsFlag
            % make sure nothing from the blacklist is included in the template training
            [allowableMotionId, allowableMotionInd] = setdiff(motionId, templateBlacklist);
            allowableMotionLocation = motionLocations(allowableMotionInd, :);
        else                
            % make sure nothing from the blacklist is included in the template training
            [allowableMotionId, allowableMotionInd] = setdiff(motionLocations(:, 1), templateBlacklist);
            allowableMotionLocation = motionLocations(allowableMotionInd, :);
        end

        
        if size(allowableMotionLocation, 1) < templateCount
            fprintf(['DataSource2: Template available less than desired for ' ...
                currMotion, ', have ' num2str(size(allowableMotionLocation, 1)) ' but want ' num2str(templateCount) '\n']);
            templateCountUse = size(allowableMotionLocation, 1);
        else
            templateCountUse = templateCount;
        end
        
        if fullSegmentsFlag
            selectIdInd = allowableMotionId(end-templateCountUse+1:end);
        else
            selectIdInd = allowableMotionLocation(end-templateCountUse+1:end, 1);
        end
        
        includeMatrix(selectIdInd, i) = 1; % NOT ID, the indice of the item itself
    end
    
    includeVector = sum(includeMatrix, 2);
    includeVector = +includeVector;
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

function templateStruct = dataSetTemplate(dataset, segmentInd, segmentLabels, segmentId, motionName, motionNameId, templateCount, frameRate, dofToUse, hmm)
    % pull out the template motion, the last 'templateCount' few

    % BW filter parameters
%     filterFreq = 0.5;
%     filterSample = 25;
%     filterOrder = 5;
    
    [motionLocations, motionId, segmentName] = motionLocator(segmentLabels, segmentId, motionName);
    
    if size(motionLocations, 1) < templateCount
        fprintf(['DataSource2: Template available less than desired for ' ...
            motionName, ', have ' num2str(size(motionLocations, 1)) ' but want ' num2str(templateCount) '\n']);
        templateCount = size(motionLocations, 1);
    end
    
    motionLocationCount = motionLocations(end-templateCount+1:end, :);
    
    templateTime = cell(1, templateCount);
    templateDataSet = cell(1, templateCount);
    segTime = zeros(1, templateCount*2);
    label = cell(1, templateCount*2);
    id = zeros(1, templateCount*2);
    
    for i = 1:templateCount
        segStart = segmentInd(motionLocationCount(i, 1), 1); % combine the "up" and "down" motion
        segEnd = segmentInd(motionLocationCount(i, end), 2);
        time = (segStart:segEnd)'/frameRate;
        
        templateData = unwrap(dataset(segStart:segEnd, dofToUse)); % strip out the frame number, and unwrap angle        
%         templateData = filter_dualpassBW(templateData, filterFreq, filterSample, filterOrder); % and filter
        
        if hmm
            time = time';
            templateData = templateData';
        end
        
        templateTime{i} = time;
        templateDataSet{i} = templateData;
        
        doubleInd = (i-1)*2+1:(i)*2;
        segTime(doubleInd) = time([1 end]);
        label{doubleInd(1)} = motionName;
        label{doubleInd(2)} = motionName;
        id(doubleInd) = [motionNameId motionNameId];

    end
    
    templateStruct.data = templateDataSet;
    templateStruct.time = templateTime;
    templateStruct.label = label;
    templateStruct.id = id;
    templateStruct.manualSegTime = segTime;
    templateStruct.count = templateCount;
end

function [segmentTime, segmentLabel, segmentCount] = dataSetSegmentation(segmentInd, segmentLabels, segmentId, motionName, segmentPassCount, frameRate)
    % pull out the template motion, the last 'templateCount' few

%     % BW filter parameters
%     filterFreq = 0.5;
%     filterSample = 25;
%     filterOrder = 5;
    
    [motionLocations, namedSegmentId, segmentName] = motionLocator(segmentLabels, segmentId, motionName);
    motionLocationInd = find(namedSegmentId(:, 1) <= segmentPassCount);
    motionLocationCount = motionLocations(motionLocationInd, :);
    segmentCount = length(motionLocationInd);
    
    segmentTime = zeros(1, segmentCount*2);
    segmentLabel = cell(1, segmentCount*2);
    
    for i = 1:length(motionLocationInd)
        segStart = segmentInd(motionLocationCount(i, 1), 1); % combine the "up" and "down" motion
        segEnd = segmentInd(motionLocationCount(i, end), 2);
%         time = (segStart:segEnd)'/frameRate;
        
        segmentTime(2*i - 1) = segStart/frameRate;
        segmentTime(2*i) = segEnd/frameRate;
        
        segmentLabel{2*i - 1} = motionName;
        segmentLabel{2*i} = motionName;
    end
end

function [motionLocations, segmentId, segmentName] = motionLocator(segmentLabels, segmentCount, desiredSegment, fullSegmentsFlag)
    % find all the occurances of the 'desiredSegment'
    motionLocations = [];
    segmentId = [];
    segmentCountExpended = []; % once the segment has been consumed, add it here
    segmentName = {};
    
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
                (strcmpi('ALLMOTIONS', desiredSegment) && ...
                (fullSegmentsFlag && ~sum(segmentCountExpended == segmentCount(i)) || ~fullSegmentsFlag)) % matching by segment id
            currSegmentCount = segmentCount(i);
            
            % if it is on fullSegment mode, and the motion is one of the
            % variants, force the overall label to be one of the accepted
            % versions
            if fullSegmentsFlag
                currLabel = fullMotionLabel(segmentLabels{i});
            else
                currLabel = segmentLabels{i};
            end
            
            if fullSegmentsFlag
                segCountInd = find(segmentCount == currSegmentCount);
            else
                segCountInd = i;
            end
            
            segmentInd = [segCountInd(1); segCountInd(end)];
            segmentCountExpended = [segmentCountExpended, currSegmentCount];
            
            if isempty(motionLocations)
                % initialize the array with 
                motionLocations = segmentInd';
                segmentId = segmentCount(i);
                segmentName{1} = currLabel;
            else
                motionLocations = [motionLocations; segmentInd'];
                segmentId = [segmentId; segmentCount(i)];
                segmentName{end+1} = currLabel;
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

function segmentData = motionAnalysis_TRO2009_segFileRead(fileToImport)
    % This function was designed to read a motionanalysis segment labelled data. It
    % was written specifically to read Dana's TRO2009 segmentation information
    
    baseFolder = 'D:\MyDocument\MotionData\TRO2009_dkulic\';
    
    if ~exist('fileToImport', 'var')
        % demo
%         fileToImport = [baseFolder, 'Jon_Rand1.trc'];
%         fileToImport = [baseFolder, 'KulicEtAlTRO_Segmentation_ManuallyLabeled.txt'];
        fileToImport = [baseFolder, 'KulicEtAlTRO_Segmentation_ManuallyLabeled_consolidated.csv'];
%         fileToImport = [baseFolder, 'KulicEtAlTRO_Segmentation_ManuallyLabeled_overlapRemoved_shortSegmentExtended.txt']; % there are some problems with the old one. this one adjusts some bounds slightly
    end
    
    delimiter = ',.'; % tab char
    
    % opening file
    fid = fopen(fileToImport,'r');
    
    % header data
    tline = fgetl(fid);
    temp = textscan(tline, '%s', 'delimiter', delimiter);
    delim = temp{1};
      
    j = 0;
    for i = 1:length(delim)
        if ~isempty(delim{i})
            j = j + 1;
            markername{j} = delim{i};
%             marker{j} = zeros(numFrame, 3);
        end
    end
    
    j = 0;
    tline = fgetl(fid);
    segmentTimeBlackListId = [];
    segmentTimeBlackListInd = [];
    
    while ischar(tline)
        temp = textscan(tline, '%s', 'delimiter', delimiter);
        delim = temp{1};
        j = j + 1;
        
        segmentId(j) = str2num(delim{1}); 
        segmentName{j} = delim{2};
        segmentTimeOrig(j, 1) = round(str2num(delim{3})/3);
        segmentTimeOrig(j, 2) = round(str2num(delim{4})/3);
        segmentTimeModified(j, 1) = round(str2num(delim{5})/3);
        segmentTimeModified(j, 2) = round(str2num(delim{6})/3);

        if segmentTimeOrig(j, 1) ~= segmentTimeModified(j, 1) || segmentTimeOrig(j, 2) ~= segmentTimeModified(j, 2)
            % if there is a difference between the original and the 
            % modified, list it for potenital blacklist
            segmentTimeBlackListId = [segmentTimeBlackListId segmentId(j)];
            segmentTimeBlackListInd = [segmentTimeBlackListInd j];
        end

        tline = fgetl(fid);
    end
    
    segmentData.header = markername;
    segmentData.segmentId = segmentId';
    segmentData.segmentName = segmentName';
    segmentData.segmentTimeOrig = segmentTimeOrig;
    segmentData.segmentTimeModified = segmentTimeModified;
    segmentData.segmentBlackList = unique(segmentTimeBlackListId);
    segmentData.segmentBlackListInd = unique(segmentTimeBlackListInd);
end

function label = fullMotionLabel(currLabel)
    switch currLabel
        case 'LKAR'
            label = 'LKE';
            
        case 'LPAR'
            label = 'LPE';
            
        case 'LPUE'
            label = 'LPE';
            
        case 'RKAR'
            label = 'RKE';
            
        case 'RPAR'
            label = 'RPE';
            
        otherwise
            label = currLabel;
    end
end