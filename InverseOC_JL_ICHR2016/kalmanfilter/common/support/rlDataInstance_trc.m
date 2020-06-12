classdef rlDataInstance_trc < rlDataInstance % < rlModel
    % hold, load, and clean input data to be passed between
    % modelinstance/featureset
    properties
        sortIndicesFlag = 0;
%         measurement = struct('name', '', ...
%             'data', [], ...
%             'use', 0); % measurement may be processed and is used in EKF
    end
    
    methods
        function obj = rlDataInstance_trc(modelInstance)
            obj = obj@rlDataInstance(modelInstance);
        end
        
        function generateVirtualMarkers(obj)
            % then use the modelInstance to generate the virtual marker data
            for i = 1:length(obj.modelInstance.jointCentreArray)
                markerData = {};
                virtualMarkerStr = obj.modelInstance.jointCentreArray{i}{1};
                for j = 2:length(obj.modelInstance.jointCentreArray{i})
                    markerStr = obj.modelInstance.jointCentreArray{i}{j};
                    
                    switch markerStr
                        case 'HARRINGTON_R'
                            [hjc, ~, ~, ~] = harrington2007prediction(obj.data, obj.modelInstance.harrington_crop1, obj.modelInstance.harrington_crop2);
                            newData = repmat(hjc', length(obj.time), 1);
                            
                        case 'HARRINGTON_L'
                            [~, hjc, ~, ~] = harrington2007prediction(obj.data, obj.modelInstance.harrington_crop1, obj.modelInstance.harrington_crop2);
                            newData = repmat(hjc', length(obj.time), 1);
                            
                        otherwise
                            newData = obj.data.(markerStr);
                    end
                    
                    markerData = [markerData newData];
                end
                
                obj.data.(virtualMarkerStr) = jointCentre(markerData);
            end
        end
        
        function loadData(obj, currFileEntry, algorithmParam)
            % use the modelInstance to load the raw data
            [obj.dt, obj.time, obj.data] = obj.modelInstance.loadData_trc(currFileEntry, algorithmParam);
            obj.generateVirtualMarkers();
            
            if 0
                % plot the full set of raw data
                vis = rlVisualizer('vis',640,480);
                
                % visualize markers
                indToPlot = 1;
                applyMarkersToVisualization(vis, [], obj, indToPlot);
                
                vis.update();
                pause(0.01);
            end
        end
        
        function plotData(obj)
            vis = rlVisualizer('vis',640,480);
            
            % visualize markers
            indToPlot = 1;
            applyMarkersToVisualization(vis, [], obj, indToPlot);
            
            vis.update();
            pause();
        end
        
        function resampleLength(obj, dt)
            obj.dt = dt;
            timeOld = obj.time;
            timeNew = timeOld(1):dt:timeOld(end);
            
            fieldNames = fieldnames(obj.data);
            for i = 1:length(fieldNames)
                dataOld = obj.data.(fieldNames{i});
                dataNew = interp1(timeOld, dataOld, timeNew);
                obj.data.(fieldNames{i}) = dataNew;
            end
            
            obj.time = timeNew;
        end
        
        function loadCrop(obj, pathToSegment)
            segmentCrop = obj.loadSegmentFile(pathToSegment);
            
            % apply cropping
            if ~isempty(segmentCrop)
                [startVal, startInd] = findClosestValue(segmentCrop.timeStart, obj.time);
                [endVal, endInd] = findClosestValue(segmentCrop.timeEnd, obj.time);
                obj.time = obj.time(startInd:endInd);
                fieldNames = fieldnames(obj.data);
                for i = 1:length(fieldNames)
                    obj.data.(fieldNames{i}) = obj.data.(fieldNames{i})(startInd:endInd, :);
                end
            end
        end

        function firstFrame = updateBaseMarker(obj, baseMarkerStr)
            if isempty(baseMarkerStr)
                firstFrame = zeros(1, 3);
                return
            end

            % subtract all the subsequent frames so that baseMarkerStr is
            % always in the same spot as the first frame
            markerOffsetVal = obj.data.(baseMarkerStr);
            firstFrame = markerOffsetVal(1, :);
%             markerOffsetVal = markerOffsetVal + repmat(firstFrame, size(obj.time, 1), 1);

            % for the purpose of the offset calculation, remove NaNs by keeping
            % it as the previous known value
            nanVals = find(isnan(markerOffsetVal(:, 1)));
            for i = 1:length(nanVals)
                lastKnownInd = nanVals(i) - 1;
                if lastKnownInd == 0
                    lastKnownValue = zeros(3, 1);
                else
                    lastKnownValue = markerOffsetVal(lastKnownInd, :);
                end
                
                markerOffsetVal(nanVals(i), :) = lastKnownValue;
            end
            
            markerOffsetVal = filter_dualpassBW(markerOffsetVal);
            
            if 0 
                figure; 
                plot(markerOffsetVal);
            end
            
            fieldNames = fieldnames(obj.data);
            for i=1:numel(fieldNames)
                obj.data.(fieldNames{i}) = obj.data.(fieldNames{i}) - markerOffsetVal;
            end
            
            for i=1:size(obj.measurement,1)
                for j=1:size(obj.measurement,2)
                    obj.measurement(i,j).mes(1:3) = obj.measurement(i,j).mes(1:3) - markerOffsetVal(i, 1:3)';
                end
            end
        end
        
        function dataProcessing(obj, algorithmParam)
            % iterate though the marker data and assign a usability label
            % add sensors
            curr_sens = SensorCore('mocapSen');
            for i = 1:length(obj.modelInstance.trcSensorDecorator)
                curr_sens.addDecorator(obj.modelInstance.trcSensorDecorator{i});
            end
            
            % run though all the combinations to which markers to add
            sensorStrToAdd = {};
%             isModelSensorToAdd = [];
            for i = 1:size(obj.modelInstance.sensorAttachmentArray, 1) % markers that is in mdl
                currSensorAttachmentArray = obj.modelInstance.sensorAttachmentArray{i};
                for j = 3:size(currSensorAttachmentArray, 2)
                    currSensorAttachment = [currSensorAttachmentArray(1) currSensorAttachmentArray(2) currSensorAttachmentArray(j)];
                    sensorStrToAdd{end+1} = currSensorAttachment{3};
%                     isModelSensorToAdd(end+1) = 1;
                end
            end
            
            if algorithmParam.addUnknownMarkersToEkf % all the other markers in the trcdata
                field_names = fieldnames(obj.data);
                for i = 1:length(field_names)
                    sensorName = field_names{i};
                    sensorNameSplit = strsplit(sensorName, '_');
                    if strcmpi(sensorNameSplit{end}, 'JC')
                        % skip this one
                        continue;
                    elseif sum(strcmpi(sensorStrToAdd, sensorName) == 1)
                        % this one has already been added
                        continue;
                    end
                    
                    sensorStrToAdd{end+1} = sensorName;
%                     isModelSensorToAdd(end+1) = 0;
                end
            end
            
            % alphabetical order
            if obj.sortIndicesFlag
                [obj.measurement_labels, sortInd] = sort(sensorStrToAdd);
            else
                obj.measurement_labels = sensorStrToAdd;
            end

%             obj.measurement_sensor = isModelSensorToAdd(sortInd);
            
            % add the marker data to the measurement array
            for i = 1:length(obj.measurement_labels)
                markerStr = obj.measurement_labels{i};
%                 isModelSensor = obj.measurement_sensor(i);
%                 fprintf('rlDataInstance: Loading sensor to dataInstance.measurement: %s \n', markerStr);
                [pos, vel, useMarker, len] = obj.checkSensorBeforeAdding(markerStr, curr_sens, algorithmParam.missingMarkerValue);
                obj.addSensorToMeas(pos, vel, useMarker, len, curr_sens, markerStr);
            end
            

%             % add the markers that are in rlModelInstance to the model and to
%             % the measurement array
%             isModelSensor = 1;
%             indMarkers = 0;
%             for i = 1:size(obj.modelInstance.sensorAttachmentArray, 1)
%                 currSensorAttachmentArray = obj.modelInstance.sensorAttachmentArray{i};
%                 for j = 3:size(currSensorAttachmentArray, 2)
%                     indMarkers = indMarkers + 1;
%                     currSensorAttachment = [currSensorAttachmentArray(1) currSensorAttachmentArray(2) currSensorAttachmentArray(j)];
%                     
%                     markerStr = currSensorAttachment{3};
%                     fprintf('rlDataInstance: Loading mdl_sensor %s...\n', markerStr);
%                     [pos, vel, useMarker, len] = obj.checkSensorBeforeAdding(markerStr, curr_sens);
%                     obj.addSensorToMeas(pos, vel, useMarker, len, curr_sens, markerStr, isModelSensor);
%                 end
%             end
%             
%             % if adding unknown markers, add the unknown markers (that does not
%             % end with '_JC' since those are synthetic) to the measurement array 
%             % so EKF can select those for markerswapping resolutions
%             if algorithmParam.addUnknownMarkersToEkf
%                 isModelSensor = 0;
%                 field_names = fieldnames(obj.data);
%                 for i = 1:length(field_names)
%                     sensorName = field_names{i};
%                     sensorNameSplit = strsplit(sensorName, '_');
%                     if strcmpi(sensorNameSplit{end}, 'JC')
%                         % skip this one
%                         continue;
%                     elseif sum(strcmpi(obj.measurement_labels, sensorName) == 1)
%                         % this one has already been added
%                         continue;
%                     end
%                     
%                     indMarkers = indMarkers + 1;
%                     
%                     markerStr = sensorName;
%                     fprintf('rlDataInstance: Loading extra_sensor %s...\n', markerStr);
%                     [pos, vel, useMarker, len] = obj.checkSensorBeforeAdding(markerStr, curr_sens);
%                     obj.addSensorToMeas(pos, vel, useMarker, len, curr_sens, markerStr, isModelSensor);
%                 end 
%             end
            
            clear curr_sens;
        end
        
%         function boolSucess = checkSensorAddMeasurement(obj, markerStr, curr_sens)
%             boolSucess = 1;
%             [pos, vel, useMarker, len] = checkSensorBeforeAdding(obj, markerStr, curr_sens);
%             
%             if useMarker == 0 || norm(pos(1, :)) == 0 
%                 % non-sensical value. don't add this sensor
%                 fprintf('checkSensorAddMeasurementInModel WARNING: %s is missing data over all frames or is missing data in the first frame and is flagged as not usable \n', markerStr);
%                 boolSucess = 0;
%             else
%                 % is okay
%             end
%             
%             isModelSensor = 1;
%             obj.addSensorToMeas(pos, vel, useMarker, len, curr_sens, markerStr, isModelSensor);
%         end
        
%         function boolSucess = checkSensorAddMeasurementInModel(obj, markerStr, curr_sens)
%             boolSucess = 1;
%             [pos, vel, useMarker, len] = checkSensorBeforeAdding(obj, markerStr, curr_sens);
%             
%             if useMarker == 0 || norm(pos(1, :)) == 0 
%                 % non-sensical value. don't add this sensor
%                 fprintf('checkSensorAddMeasurementInModel WARNING: %s is missing data over all frames or is missing data in the first frame and is flagged as not usable \n', markerStr);
%                 boolSucess = 0;
%             else
%                 % is okay
%             end
%             
%             isModelSensor = 1;
%             obj.addSensorToMeas(pos, vel, useMarker, len, curr_sens, markerStr, isModelSensor);
%         end
%         
%         function boolSucess = checkSensorAddMeasurementNotInModel(obj, markerStr, curr_sens)
%             boolSucess = 1;
%             [pos, vel, useMarker, len] = checkSensorBeforeAdding(obj, markerStr, curr_sens);
%             
%             if useMarker == 0
%                 % non-sensical value. don't add this sensor
%                 fprintf('checkSensorAddMeasurementNotInModel WARNING: %s is missing data over all frames and is flagged as not usable \n', markerStr);
%                 boolSucess = 0;
%             else
%                 % is okay
%             end
%             
%             isModelSensor = 0;
%             obj.addSensorToMeas(pos, vel, useMarker, len, curr_sens, markerStr, isModelSensor);
%         end
        
        function [pos, vel, useMarker, len] = checkSensorBeforeAdding(obj, markerStr, curr_sens, missingMarkerValue)
            markerData = obj.data.(markerStr);
            len = length(obj.time);
            
            pos = markerData;
            vel = [reshape(filter_dualpassBW(diff(pos)),[],3)/obj.dt;[0 0 0]];
            
            switch missingMarkerValue
                case 'same'
                    %Create a measurement that keeps position constant       
                    pos_mes_this_marker = pos;
                    missing_markers = ismember(pos_mes_this_marker,[0 0 0],'rows');
                    %Note this actually finds the index right before missing
                    missing_starts = strfind(missing_markers',[0 1]);
                    missing_ends = strfind(missing_markers',[1 0]);
%                     if(numel(missing_starts) > numel(missing_ends))
                     if(missing_markers(end))
                        missing_ends = [missing_ends size(pos_mes_this_marker,1)];
                    end
                    if(missing_markers(1))
                        missing_starts = [1 missing_starts];
                    end
                    for(k=1:numel(missing_starts))
                        pos_mes_this_marker(missing_starts(k):missing_ends(k),:) = ...
                            repmat(pos_mes_this_marker(missing_starts(k),:),missing_ends(k)-missing_starts(k)+1,1);
                    end
                    
                    pos = pos_mes_this_marker;
                    
                    
                otherwise
                    missing_markers = pos == 0;
                    
                    switch curr_sens.binary_type
                        case 1
                            
                        case 17 % sensors is pos/vel markers
                            % apply dilation to the marker
                            missing_markers = imdilate(missing_markers, ones(25));
                    end
                    
                    if norm(pos(1, :)) > 0
                        missing_markers(1, :) = 0; % if the first frame has value, don't blank it out
                    end
                    pos(missing_markers) = missingMarkerValue;
            end
            
            useMarker = curr_sens.binary_type;
        end
        
        function addSensorToMeas(obj, pos, vel, useMarker, len, curr_sens, markerStr)
            curr_mes = SensorMeasurement(len, 1);
            [curr_mes.type] = deal(useMarker);
            [curr_mes.size] = deal(numel(curr_sens.measurement));
            
            switch useMarker
                case 1
                    data = [pos];
                    
                case 17
                    data = [pos vel];
            end
            
            curr_mes.setMesArray(data);
            
            if isempty(obj.measurement)
                % first entry in the measurement arrays
                obj.measurement = curr_mes;
%                 obj.measurement_labels{1} = markerStr;
%                 obj.measurement_sensor(end+1) = isModelSensor;
            else
                % not the first entry in the measurement arrays
                obj.measurement = [obj.measurement curr_mes];
%                 obj.measurement_labels{end+1} = markerStr;
%                 obj.measurement_sensor(end+1) = isModelSensor;
            end
        end
        
        function ind = findSensorFrame(obj, markerStr)
             ind = find(ismember(obj.measurement_labels, markerStr) == 1);
        end
        
        function currMeasurement = findMeasurement(obj, markerStr)
            if ~iscell(markerStr)
                markerStrCell{1} = markerStr;
            else
                markerStrCell = markerStr;
            end
            
            for i = 1:length(markerStrCell)
                indFrame = obj.findSensorFrame(markerStrCell{i});
                currMeasurement(:, i) = obj.measurement(:, indFrame);
            end
        end
        
        function boolUse = checkUseMarker(obj, markerStr)
            % check the use array to see if we should use a given marker or not
            indFrame = obj.findSensorFrame(markerStr);
            
            if ~isempty(indFrame)
                boolUse = obj.measurement(indFrame).use;
            else
                fprintf('WARNING: checkUseMarker did not find frame: %s\n', markerStr);
                boolUse = -1;
            end
        end
        
%         function dataArray = outputTrcMask(obj, modelInstance, featureSet, dataInstance, matchMatrix, output_matching)
% %              ekf_eventhandler = featureSet.ekf_eventhandler;
%             ekf_match = featureSet.ekf_match;
% %             ekf_labels = obj.measurement_labels;
% %             all_labels = dataInstance.measurement_labels;
%             
%             % ekf_match: len is mdl.sensors, vals is indices of meas.allSens
%             % matchMatrix: len is mdl.sensors, vals is indices of meas.allSen
%             % outputMatch: len is mdl.allSens, vals is indices of mdl.sensor
%             dataArray = zeros(length(obj.time), length(output_matching));
%             for i = 1:length(output_matching)
%                 currMes = output_matching(i);
%                 if currMes ~= 0
%                     ekfErrorEntries = ekf_match(:, currMes) == matchMatrix(i);
%                     dataArray(ekfErrorEntries, i) = 1; % maintained in the same order as ekf sensor (orig)
%                 end
%             end
%         end
    end
end