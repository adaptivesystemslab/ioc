classdef rlDataInstance_imu < rlDataInstance % < rlModel
    % hold, load, and clean input data to be passed between
    % modelinstance/featureset
    properties

    end
    
    methods
        function obj = rlDataInstance_imu(modelInstance)
            obj = obj@rlDataInstance(modelInstance);
        end
        
        function loadData(obj, currFileEntry, algorithmParam)
            [obj.data] = obj.modelInstance.loadData_imu(currFileEntry, algorithmParam);
            
            % some imus have repeated timestamps. remove these entries
            for i = 1:length(obj.data)
                [timeOld, uniqueInd] = unique(obj.data(i).time);
                gyroOld = obj.data(i).gyroscopeCalibrated(uniqueInd, :);
                accelOld = obj.data(i).accelerometerCalibrated(uniqueInd, :);
                
                % apply gyro bias
                gyroOffset = mean(gyroOld(1:200, :), 1);

                timeNew = timeOld;
                gyroNew = gyroOld - gyroOffset;
                accelNew = accelOld;

                obj.data(i).time = timeNew;
                obj.data(i).gyroscopeCalibrated = gyroNew;
                obj.data(i).accelerometerCalibrated = accelNew;
                obj.data(i).dt = mean(diff(obj.data(i).time));
            end
            
            if 0
               figure; 
               for i = 1:5
                   ax(i) = subplot(2, 5, i);
                   plot(obj.data(i).accelerometerCalibrated);
                   title(['Accel_' obj.data(i).name]);
                   ax(i+5) = subplot(2, 5, i+5);
                   plot(obj.data(i).gyroscopeCalibrated);
                   title(['Gyro_' obj.data(i).name]);
               end
               
               linkaxes(ax(1:5));
               linkaxes(ax(6:10));
            end
        end
        
        function resampleLength(obj, dt)
            % unify the time array  
            startTime = 0;
            endTime = inf;
            for i=1:numel(obj.data)
                if obj.data(i).time(1) > startTime
                    startTime = obj.data(i).time(1);
                end
                
                if obj.data(i).time(end) < endTime
                    endTime = obj.data(i).time(end);
                end
            end
            timeNew = startTime:dt:endTime;
            
            for i = 1:length(obj.data)
                timeOld = obj.data(i).time;
                dataOld = obj.data(i).gyroscopeCalibrated;
                dataNew = interp1(timeOld, dataOld, timeNew);
                obj.data(i).gyroscopeCalibrated = dataNew;
                
                timeOld = obj.data(i).time;
                dataOld = obj.data(i).accelerometerCalibrated;
                dataNew = interp1(timeOld, dataOld, timeNew);
                obj.data(i).accelerometerCalibrated = dataNew;
                
                obj.data(i).time = timeNew;
            end
            
            obj.dt = dt;
            obj.time = timeNew;
        end
        
        function loadCrop(obj, pathToSegment)
            segmentCrop = obj.loadSegmentFile(pathToSegment);
            
            % apply cropping
            if ~isempty(segmentCrop)
                for i = 1:length(obj.data)
                    [startVal, startInd] = findClosestValue(segmentCrop.timeStart, obj.data(i).time);
                    [endVal, endInd] = findClosestValue(segmentCrop.timeEnd, obj.data(i).time);
                
                    obj.data(i).time =  obj.data(i).time(startInd:endInd);
                    obj.data(i).gyroscopeCalibrated = obj.data(i).gyroscopeCalibrated(startInd:endInd, :);
                    obj.data(i).accelerometerCalibrated = obj.data(i).accelerometerCalibrated(startInd:endInd, :);
                end
            end
        end
        
        function dataProcessing(obj, algorithmParam)
            % iterate though the marker data and assign a usability label
            % add sensors
            curr_sens = SensorCore('imuSen');
            for i = 1:length(obj.modelInstance.imuSensorDecorator)
                curr_sens.addDecorator(obj.modelInstance.imuSensorDecorator{i});
            end
            
            if obj.modelInstance.imuSensorYaw
                yaw_sens = SensorCore('yawSen');
                yaw_sens.addDecorator('yaw');
            else
                yaw_sens = [];
            end
            
%             [obj.measurement, ~, obj.measurement_labels, obj.time] = loadDataFromImu(obj.data, curr_sens, yaw_sens);
            
%             indMarkers = 0;
%             for i = 1:size(obj.modelInstance.sensorImuAttachmentArray, 1)
%                 currSensorAttachmentArray = obj.modelInstance.sensorImuAttachmentArray{i};
%             end
            
            [obj.measurement, obj.data, obj.measurement_labels, obj.time] = loadDataFromImu(obj.data, curr_sens, yaw_sens);
            
            obj.dt = mean([obj.data.dt]);
            clear curr_sens yaw_sens;
        end
    end
end