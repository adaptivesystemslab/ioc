classdef rlDataInstance < handle
    properties
        modelInstance;
        
        dt;
        time;
        
        data; % raw stuff loaded from file
        measurement; % modified data to fit EKF usages
        measurement_labels;
        measurement_sensor; % 1 if part of sensor, 0 if not

        segments;
    end
    
    methods
        function obj = rlDataInstance(modelInstance)
            obj.modelInstance = modelInstance;
        end
        
        function loadSegment(obj, pathToSegment)
            obj.segments = obj.loadSegmentFile(pathToSegment);
        end
        
        function segment = loadSegmentFile(obj, pathToSegment)
            pathToRawHeader = strrep(pathToSegment, '.data', '.header');
            pathToRawData = pathToSegment;
            segmentRaw = segmentationDataHandle_crop(pathToRawHeader, pathToRawData);
            segmentRaw.load();
            
            segment = [];
            for i = 1:length(segmentRaw.timeStart) % change it to array of structs for easy indexing
                segment(i).timeStart = segmentRaw.timeStart(i);
                segment(i).timeEnd = segmentRaw.timeEnd(i);
            end
        end
        
        function updateData(obj, time, data)
            obj.time = time;
            obj.dt = mean(diff(time));
            obj.data = data;
        end
        
        function updateMeasurement(obj, time, measurement, labels)
            if ~isempty(time)
                obj.time = time;
                obj.dt = mean(diff(time));
            end
            
            % check to see if the new measurements contain the same type of
            % data as the old measurements (but potentially in a different
            % order)
            sameFlag = 1;
            for i = 1:length(labels)
                labelIndCurr = find(ismember(obj.measurement_labels, labels{i}));
                if isempty(labelIndCurr)
                   sameFlag = 0;
                   break;
                else
                    labelInd(i) = labelIndCurr;
                end
            end
            
            if sameFlag
                % same entries, but potentially could be out of order
                newMes = SensorMeasurement(size(measurement, 1), size(measurement, 2));
                for i = 1:length(labels)
                    labelInd = find(ismember(obj.measurement_labels, labels{i}));
                    newMes(:, labelInd) = measurement(:, i);
                end
            else
                % different data, just replace it
                obj.measurement = measurement;
                obj.measurement_labels = labels;
            end
        end
        
        function matchMatrix = outputMatchMatrix(obj, mdl)
            % find the entries that match the name. the first col is for the
            % mdl, the second col is for meas
            matchMatrix = zeros(length(obj.measurement_labels), 2);
            matchMatrix(1:length(mdl.sensors)) = 1:length(mdl.sensors);
            
            % move all the entries that is actively in the model to the front
            sensorNames = {mdl.sensors.name};
            for i = 1:length(sensorNames)
                indx = find(ismember(obj.measurement_labels, sensorNames{i}), 1);
                if ~isempty(indx) && matchMatrix(i, 2) == 0
                    matchMatrix(i, 2) = indx;
                else
                    lala = 1;
                end
            end
            
            % then all the rest
            leftoverInds = find(ismember(1:length(obj.measurement_labels), matchMatrix(1:length(sensorNames), 2)') == 0);
            matchMatrix(length(sensorNames)+1:length(obj.measurement_labels), 2) = leftoverInds;
        end
        
        function [time, data, indUse, label, matchMatrix] = outputMeasurementForEkf(obj, mdl, indUse)
            if ~exist('indUse', 'var')
                indUse = 1:length(obj.time);
            end
            
            % based on the model, rearrange the measurement output so that the
            % first 'n' entries matches the order of the mdl, and the last 'm'
            % entries are the un-modelled data (if approprate)
            matchMatrix = obj.outputMatchMatrix(mdl);

            time = obj.time(indUse);
            data = obj.measurement(indUse, :);
            label = obj.measurement_labels;
        end
        
        function [time, data, indUse, label, matchMatrix] = repeatMeasurementForEkf(obj, mdl, indFrameToRepeat, numRepeat)
            % based on the model, rearrange the measurement output so that the
            % first 'n' entries matches the order of the mdl, and the last 'm'
            % entries are the un-modelled data (if approprate)
            matchMatrix = obj.outputMatchMatrix(mdl);
%             matchMatrix = [1:length(obj.measurement_labels); 1:length(obj.measurement_labels)]';
            
            time = (0:numRepeat)*obj.dt;

            % duplicate the target array
            data = obj.rowToRepeat(obj.measurement(indFrameToRepeat, :), numRepeat);
            indUse = indFrameToRepeat*ones(size(time));
            
            % rearrange the target array
%             data = data_tmp(:, matchMatrix(:, 2));
%             label = obj.measurement_labels(matchMatrix(:, 2));

            label = obj.measurement_labels;
        end
        
        function data = rowToRepeat(obj, dataRow, numRepeat)
            data_tmp = SensorMeasurement(numRepeat, size(dataRow, 2));
            indUse = zeros(1, numRepeat);
            
            for i = 1:numRepeat
                data_tmp(i, :) = dataRow;
            end
            
            data = data_tmp;
        end
        
        function matchMatrix = outputMatchMatrixCheckDataInstance(obj, mdl, modelInstance)
             % find the entries that match the name. the first col is for the
            % mdl, the second col is for meas
            matchMatrix = zeros(length(obj.measurement_labels), 2);
%             matchMatrix(1:length(mdl.sensors)) = 1:length(mdl.sensors);
            
            % move all the entries that is actively in the model to the front
            sensorNames = {mdl.sensors.name};
            addInd = 0;
            for i = 1:length(sensorNames)
                indx = find(ismember(obj.measurement_labels, sensorNames{i}), 1);
                if ~isempty(indx) && matchMatrix(i, 2) == 0
                    % now check for sensor data integrety
                    lengthUseInd = modelInstance.lengthUseInd;
%                     sensorMarker = obj.measurement(lengthUseInd, indx).getMesArray;
                    sensorMarker = obj.data.(sensorNames{i});
                    bool = checkSensorPlacement([], sensorMarker, [], [], [], lengthUseInd, '', '');
                    
                    if bool
                        addInd = addInd + 1;
                        matchMatrix(addInd, 1) = i;
                        matchMatrix(addInd, 2) = indx;
                    else
                        lala = 1;
                    end
                else
                    lala = 1;
                end
            end
            
            % then all the rest
            leftoverInds = find(ismember(1:length(obj.measurement_labels), matchMatrix(1:addInd, 2)') == 0);
            matchMatrix(addInd+1:length(obj.measurement_labels), 2) = leftoverInds;
        end
        
        function removeMeasurements(obj, sensorList)
            if ~isempty(sensorList)
                labelKeepInds = zeros(1, length(obj.measurement_labels));
                for i = 1:length(obj.measurement_labels)
                    labelIndCurr = find(ismember(sensorList, obj.measurement_labels{i}));
                    if isempty(labelIndCurr)
                        labelKeepInds(i) = 0;
                    else
                        labelKeepInds(i) = i;
                    end
                end
                
                labelKeepInds = labelKeepInds(labelKeepInds > 0);
                
                obj.measurement = obj.measurement(:, labelKeepInds);
                obj.measurement_labels = obj.measurement_labels(labelKeepInds);
            end
        end
    end
end