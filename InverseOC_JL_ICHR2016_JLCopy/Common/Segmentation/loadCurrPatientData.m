function [subjectDataTemp, segmentCropOutput, segmentManualOutput, ekfOutput, imuOutput, forceOutput, primitiveName, useInTraining, mocapOutput] = ...
    loadCurrPatientData(obj, currSubjInfo)
    % load patient data, and separate out the segmentation and joint angles
    % variables
    
    % set up blank variables
    subjectDataTemp = [];
    segmentCropLoaded = [];
    segmentManualLoaded = [];
    ekfLoaded = [];
    imuLoaded = [];
    forceLoaded = [];

    % load the full package first
    try
        
        fprintf('    Loading patient data: %s\n', currSubjInfo.filePath);
%         switch currSubjInfo.datasetName
%             case {'squats_tuat_2011', 'squats_tuat_2015', 'doppel', 'taiso_ut_2009'}
                subjectDataTemp = exerciseDataHandle_generalized(currSubjInfo);
                subjectDataTemp.load;
                
%             otherwise % APARS file format
%                 subjectDataTemp = exerciseDataHandle(currSubjInfo);
%                 subjectDataTemp.load;
%         end

    catch err
        if strcmpi(obj.settings.mode, 'Training')
            % if we're in training mode, skip this instance
            fprintf('Training (loading) - ignoring load error: %s\n', err.message);
            return
        else
            % if we're in testing mode, break out of the function
            rethrow(err);
        end
    end
    
    % loading segmentation information
    try
        fprintf('    Loading segment: %s\n', currSubjInfo.filepathDataSegManual);

        switch currSubjInfo.datasetName
            case {'squats_tuat_2011', 'squats_tuat_2015', 'doppel', 'taiso_ut_2009'}
                segmentCropLoaded = subjectDataTemp.dataSegCrop;
                segmentManualLoaded = subjectDataTemp.dataSegManual; % well...the class is too rigid for the purposes of this fct
                
            otherwise % APARS file format
                segmentCropLoaded = subjectDataTemp.dataSegCrop;
                segmentManualLoaded = subjectDataTemp.dataSegManual.convertToStruct; % well...the class is too rigid for the purposes of this fct
        end
        
        if isempty(segmentManualLoaded) || length(segmentManualLoaded.timeStart) < obj.segmentLengthThreshold % if there is not enough segments available, just skip it
            fprintf('Insufficient segment samples: %s', currSubjInfo.filepathDataSegManual);
        end

    catch err
        if strcmpi(obj.settings.mode, 'Training')
            % if we're in training mode, skip this instance
            fprintf('Training (seg) - ignoring load error: %s\n', err.message);
            return
        else
            % if we're in testing mode, break out of the function
            rethrow(err);
        end
    end


    % loading EKF and joint angles
    try
        fprintf('    Loading data streams: %s\n', currSubjInfo.filepathDataEkf);

        switch currSubjInfo.datasetName
            case 'wojtusch2015_ichr'
                mocapLoaded = [];
                imuLoaded = [];
                [ekfLoaded, forceLoaded] = loadWojtusch2015_ICHR_mat(currSubjInfo.ekfPath);
                
            case 'doppel'
                mocapLoaded = subjectDataTemp.dataMocap;
                imuLoaded = [];
                ekfLoaded = subjectDataTemp.dataEkf;
                forceLoaded = [];
                
            case {'squats_tuat_2011', 'squats_tuat_2015', 'taiso_ut_2009'}
                mocapLoaded = subjectDataTemp.dataMocap;
                imuLoaded = [];
                ekfLoaded = subjectDataTemp.dataEkf;
                forceLoaded = subjectDataTemp.dataForces;
                
            otherwise
                % normal loading
                mocapLoaded = [];
                imuLoaded = subjectDataTemp.dataImu(1);
                ekfLoaded = subjectDataTemp.dataEkf;
                forceLoaded = [];
        end

        if length(ekfLoaded.time) < obj.ekfLengthThreshold
            if strcmpi(obj.settings.mode, 'Training')
                % if we're in training mode, skip this instance
                return
            else
                % if we're in testing mode, break out of the function
                throw('Improperly shaped EKF time data');
            end
        end

    catch err
        if strcmpi(obj.settings.mode, 'Training')
            % if we're in training mode, skip this instance and
            % try loading the next training template
            fprintf('Training (ekf/imu) - ignoring load error: %s', err.message);
            return
        else
            % if we're in testing mode, break out of the
            % function and continue to the next test
            rethrow(err);
        end
    end
    
    % examine the 'segmentManualLoaded'. if loaded and has multiple
    % subjects attributed to it, then split them into smaller pieces
    % according to some metric. 
    
    % regardless, the first item in the cell will be retained as the
    % full length
    segmentCropOutput{1} = segmentCropLoaded;
    segmentManualOutput{1} = segmentManualLoaded;
    ekfOutput{1} = ekfLoaded;
    imuOutput{1} = imuLoaded;
    forceOutput{1} = forceLoaded;
    mocapOutput{1} = mocapLoaded;
    primitiveName{1} = subjectDataTemp.exerciseType;
    useInTraining{1} = 1;
    
    if exist('segmentManualLoaded', 'var')
        switch class(segmentManualLoaded)
            case 'struct'
                if ~isfield(segmentManualLoaded, 'primitiveId') || isempty(segmentManualLoaded.primitiveId)
                    % there is no 'primitiveID' defined, so no splitting is needed
                    return
                end
            case 'segmentationDataHandle'
                if isempty(segmentManualLoaded.primitiveId)
                    % there is no 'primitiveID' defined, so no splitting is needed
                    return
                end
        end
        
        uniquePrimitives = unique(segmentManualLoaded.primitiveId);
        if length(uniquePrimitives) == 1
            % the segments are all the same type of motion, so we don't
            % need to proceed with the splitting
            
            % Doppel case with one bodyID in segments but multiple in input
            % data
            count = 1;
            if(strcmp(currSubjInfo.datasetName, 'doppel'))
                ID = uniquePrimitives;
                
                index = find(ekfLoaded.bodyID == ID);
%                 primitiveName{1} = ['Body ID ' num2str(ID)]; %TODO will need to populate this with some sort of identifier. subject string probably
                         
                %Separate by Body ID first, then remove duplicates
                dataNames = fieldnames(ekfLoaded);
                for name = 1:length(dataNames)
                    % Save all loaded fields for this body ID
                    tempData.(dataNames{name}) = ekfLoaded.(dataNames{name})(index,:);
                end

                time = tempData.time;
                uniqueTime = unique(time);
                
%                 if(length(time) ~= length(uniqueTime))
                    indices = [];
                    for k = 1:length(uniqueTime)
                        tempind = find(tempData.time == uniqueTime(k), 1); %take only the first one
                        indices = [indices tempind];
                    end
%                 else
%                     indices = index;
%                 end
                
                for name = 1:length(dataNames)
                    % Save all loaded fields for this body ID
                    curData.(dataNames{name}) = tempData.(dataNames{name})(indices,:);
                end
                
                
                

                %Compute and save angle velocities
                [len wid] = size(diff(curData.Q));
                timediff = (repmat(diff(curData.time),1,wid));
                % Derivative of angle data
                curData.dQ = [zeros(1, size(curData.Q, 2)); ...
                    diff(curData.Q)./(timediff)];
                
                
                % Get segment indices for current bodyId
                inds = find(segmentManualLoaded.primitiveId == ID);
                
                
                
                
                
                curSeg.timeStart = segmentManualLoaded.timeStart(inds);
                curSeg.timeEnd = segmentManualLoaded.timeEnd(inds);
                curSeg.primitiveId = segmentManualLoaded.primitiveId(inds);
                %curSeg.segmentId = segmentManualLoaded.segmentId(inds);
                curSeg.segmentName = segmentManualLoaded.segmentName(inds);
                curSeg.segmentCount = (1:length(curSeg.timeStart))';
                curSeg.use = segmentManualLoaded.use(inds);
                
                  % Do we want to split this ID further at large time
                    % discontinuities?
                    if(obj.splitAtTimeGaps)
                        [splitData, splitSegs] = SplitAtTimeGaps(obj, curData, curSeg);
                        
                        
                        
                        % Save each of the split data individually
                        for j = 1:length(splitData)
                            Data = splitData{j};
                            Seg = splitSegs{j};
                            
                            if(obj.noJumpRope > 0) %Remove all skipping sections                                
                                try
                                if(strcmp(Seg.segmentName(1), 'JumpRope'))
                                    continue;
                                end
                                catch
                                end
                            elseif(obj.noJumpRope < 0) %Remove all non-skipping sections
                                try
                                if(~(strcmp(Seg.segmentName(1), 'JumpRope')))
                                    continue;
                                end
                                catch
                                end
                                
                            end
                            
                            segmentCropOutput{count} = [];
                            segmentManualOutput{count} = Seg;
                            ekfOutput{count} = Data;
                            imuOutput{count} = [];
                            forceOutput{count} = [];
                            
                            useInTraining{count} = 1;
                            primitiveName{count} = ['Section Number ' num2str(count) ' Body ID ' num2str(ID) '-' Seg.segmentName{1}];
                            
                            count = count+1;
                        end
                    else
%                         Save extracted values (don't what to keep full
%                         length file in index 1)
                            segmentCropOutput{count} = [];
                            segmentManualOutput{count} = curSeg;
                            ekfOutput{count} = curData;
                            imuOutput{count} = [];
                            forceOutput{count} = [];
                            
                            useInTraining{count} = 1;
                            primitiveName{count} = ['Body ID ' num2str(ID)];
                            
                            count = count+1;
                            
                    end
                
                % Save extracted values (don't what to keep full
                % length file in index 1)
%                 segmentCropOutput{1} = [];
%                 segmentManualOutput{1} = curSeg;
%                 ekfOutput{1} = curData;
%                 imuOutput{1} = [];
%                 forceOutput{1} = [];
%                 
%                 useInTraining{1} = 1;
                
                h = figure;
                hold on, plot(curData.time, curData.jp(:,2)), plotBoxes(h, curSeg.timeStart, curSeg.timeEnd);
                 title(['Body ID: ' num2str(ID) ': Input Segments and Y-position of ' obj.jpNames{1}])
                ylabel('Joint Position [m]');
                xlabel('Raw time');
                
                fileName = ['Loaded Segments - BodyID', num2str(ID), '-Subject', num2str(currSubjInfo.subjectNumber)];
                fileName = [fileName, '.fig'];
                path = obj.figureExportPath;
                toSave = fullfile(path, fileName);
                
                e = exist(path, 'dir');
                if(e)
                    if exist(toSave, 'file')
                       
                    end
                    saveas(h, toSave);
                else
                    mkdir(path)
                    saveas(h,toSave);
                end
                
                close(h);
                return;
            end
            

            return
        end
        
        switch currSubjInfo.datasetName
            case {'taiso_ut_2009'}
                useInTraining{1} = 0; % don't use the full length data for training if we're going to include the individuals
                  
                for ind_uniquePrimitives = 1:length(uniquePrimitives)
                    % separate the segment file according to the different
                    % primitiveids
                    currSegmentManual = segmentManualLoaded.extractSubsegmentation(uniquePrimitives(ind_uniquePrimitives));
                    
                    % depending on the mode we are in, we may choose to
                    % keep different segments
                    switch obj.settings.mode
                        case 'Training'
                            currSegmentManual.use_applyTrain;
                            
                        case 'Testing'
                            currSegmentManual.use_applyTest;
                    end
                    
                    currSegmentManual = currSegmentManual.dropUseZeroSegments;
                    
                    segmentCropOutput{ind_uniquePrimitives+1} = segmentCropLoaded;
                    segmentManualOutput{ind_uniquePrimitives+1} = currSegmentManual;
                    ekfOutput{ind_uniquePrimitives+1} = ekfLoaded;
                    imuOutput{ind_uniquePrimitives+1} = imuLoaded;
                    forceOutput{ind_uniquePrimitives+1} = forceLoaded;
                    mocapOutput{ind_uniquePrimitives+1} = mocapLoaded;
                        
                    % also note some identifying feature for this set
                    strArr = strsplit(currSegmentManual.segmentName{1}, '_');
                    primitiveName{ind_uniquePrimitives+1} = strArr{1};
                    
                    useInTraining{ind_uniquePrimitives+1} = 1;
                end
                
            case {'doppel'}
               bodyIDs = uniquePrimitives;
               count = 1;  
               
                % Split input data and segments by body ID identifier
                for ID = 1:length(bodyIDs)
                    % Get indices for each body ID
                    index = find(ekfLoaded.bodyID == bodyIDs(ID));
 %TODO will need to populate this with some sort of identifier. subject string probably
                    
                    %Separate by Body ID first, then remove duplicates
                    dataNames = fieldnames(ekfLoaded);
                    for name = 1:length(dataNames)
                        % Save all loaded fields for this body ID
                        tempData.(dataNames{name}) = ekfLoaded.(dataNames{name})(index,:);
                    end
                    
                    time = tempData.time;
                    uniqueTime = unique(time);
                    
                    indices = [];
                    for k = 1:length(uniqueTime)
                        tempind = find(tempData.time == uniqueTime(k), 1); %take only the first one
                        indices = [indices tempind];
                    end
                    
                    
                    for name = 1:length(dataNames)
                        % Save all loaded fields for this body ID
                        curData.(dataNames{name}) = tempData.(dataNames{name})(indices,:);
                    end
                    
                    
                    %Compute and save angle velocities
                    [len wid] = size(diff(curData.Q));
                    timediff = (repmat(diff(curData.time),1,wid));
                    % Derivative of angle data
                    curData.dQ = [zeros(1, size(curData.Q, 2)); ...
                        diff(curData.Q)./(timediff)];
                    
                    
                    % Get segment indices for current bodyId
                    inds = find(segmentManualLoaded.primitiveId == bodyIDs(ID));
                    
                    curSeg.timeStart = segmentManualLoaded.timeStart(inds);
                    curSeg.timeEnd = segmentManualLoaded.timeEnd(inds);
                    curSeg.primitiveId = segmentManualLoaded.primitiveId(inds);
                    %curSeg.segmentId = segmentManualLoaded.segmentId(inds);
                    curSeg.segmentName = segmentManualLoaded.segmentName(inds);
                    curSeg.segmentCount = (1:length(curSeg.timeStart))';
                    curSeg.use = segmentManualLoaded.use(inds);
                    
                    % Do we want to split this ID further at large time
                    % discontinuities?
                    if(obj.splitAtTimeGaps)
                        [splitData, splitSegs] = SplitAtTimeGaps(obj, curData, curSeg);
                        
                        
                        
                        % Save each of the split data individually
                        for j = 1:length(splitData)
                            Data = splitData{j};
                            Seg = splitSegs{j};
                            
                            if(obj.noJumpRope > 0) %Remove all skipping sections                                
                                try
                                    if(strcmp(Seg.segmentName(1), 'JumpRope'))
                                        continue;
                                    end
                                catch
                                end
                            elseif(obj.noJumpRope < 0) %Remove all non-skipping sections
                                try
                                    if(~(strcmp(Seg.segmentName(1), 'JumpRope')))
                                        continue;
                                    end
                                catch
                                end
                                
                            end
                            
                        
                             %Save the section
                            segmentCropOutput{count} = [];
                            segmentManualOutput{count} = Seg;
                            ekfOutput{count} = Data;
                            imuOutput{count} = [];
                            forceOutput{count} = [];
                            
                            useInTraining{count} = 1;
                            try
                            primitiveName{count} = ['Section Number ' num2str(count) ' Body ID ' num2str(bodyIDs(ID)) '-' Seg.segmentName{1}];
                            catch
                                primitiveName{count} = ['Section Number ' num2str(count) ' Body ID ' num2str(bodyIDs(ID)) '-No Segments'];
                            end
                            
                            count = count+1;
                        end
                    else
%                         Save extracted values (don't what to keep full
%                         length file in index 1)
                            segmentCropOutput{count} = [];
                            segmentManualOutput{count} = curSeg;
                            ekfOutput{count} = curData;
                            imuOutput{count} = [];
                            forceOutput{count} = [];
                            
                            useInTraining{count} = 1;
                            primitiveName{count} = ['Body ID ' num2str(bodyIDs(ID))];
                            
                            count = count+1;
                            
                    end
                    

                    h = figure;
                    hold on, plot(curData.time, curData.jp(:,2)), plotBoxes(h, curSeg.timeStart, curSeg.timeEnd);
                    title(['Body ID: ' num2str(bodyIDs(ID)) ': Input Segments and Y-position of ' obj.jpNames{1}])
                    ylabel('Joint Position [m]');
                    xlabel('Raw time');
                    
                    fileName = ['Loaded Segments - BodyID', num2str(bodyIDs(ID)), '-Subject', num2str(currSubjInfo.subjectNumber)];
                    fileName = [fileName, '.fig'];
                path = obj.figureExportPath;
                toSave = fullfile(path, fileName);
                
                e = exist(path, 'dir');
                if(e)
                    if exist(toSave, 'file')
                        
                    end
                    saveas(h, toSave);
                else
                    mkdir(path)
                    saveas(h,toSave);
                end
                    
                    close(h)
                    
                    
                end

                
        end
    end
end