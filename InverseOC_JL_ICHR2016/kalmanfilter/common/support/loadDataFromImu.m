function [mes, imus, labels, timeOut] = loadDataFromImu(imus, curr_sens, yaw_sens)
    % find the shortest length of time array accounting for the dataInstance
    % mocap data
%     dt = mean([imus.dt]);
    
%     mesStruct = struct('name', '', ...
%         'data', [], ...
%         'use', 0); % measurement may be processed and is used in EKF
    
%     if exist('dataInstance', 'var')
% %         timeMocap = dataInstance.time;
% %         startTime = timeMocap(1);
% %         endTime = timeMocap(end);
%         startTime = 0;
%         endTime = inf;
%         
%         for i=1:numel(imus)
%             if imus(i).time(1) > startTime
%                 startTime = imus(i).time(1);
%             end
%             
%             if imus(i).time(end) < endTime
%                 endTime = imus(i).time(end);
%             end
%         end
%         
%         timeNew = startTime:dt:endTime;
%         
%         % now reinterpolate the imus to match this new array
%         for i = 1:length(imus)
%             [timeOld, uniqueInd] = unique(imus(i).time);
%             gyroOld = imus(i).gyroscopeCalibrated(uniqueInd, :);
%             accelOld = imus(i).accelerometerCalibrated(uniqueInd, :);
%             
%             gyroNew = interp1(timeOld, gyroOld, timeNew);
%             accelNew = interp1(timeOld, accelOld, timeNew);
%             
%             imus(i).time = timeNew;
%             imus(i).gyroscopeCalibrated = gyroNew;
%             imus(i).accelerometerCalibrated = accelNew;
%         end
%     end

    labels = {};
    
    %Crop IMUs
    min_samples = inf;
    for i=1:numel(imus)
        if(numel(imus(i).time) < min_samples)
            min_samples = numel(imus(i).time);
            timeOut = imus(i).time;
        end
    end
    %Build sensor measurement array we can just toss into EKF
    mes = SensorMeasurement(min_samples,numel(imus));
    %Set type as Gyro+Accel, we assume thats what the model already
    %thinks it is
    [mes.type] = deal(curr_sens.binary_type);
    [mes.size] = deal(numel(curr_sens.measurement));
    for i=1:numel(imus)
        data = [imus(i).gyroscopeCalibrated(1:min_samples,:)...
            imus(i).accelerometerCalibrated(1:min_samples,:)];
        
        mes(:,i).setMesArray(data);
        
        labels{end+1} = imus(i).name;
        
%         newMesStruct.name = imus(i).name;
%         newMesStruct.data = data;
%         newMesStruct.use = 1;
%         mesStruct(i) = newMesStruct;
    end

    %Add Yaw sensor measurements if there are Yaw sensors, Note
    %this stuff assumes yaw sensors are attached to the model after
    %actual IMUs
%     for i=1:numel(mdl.sensors)
    if ~isempty(yaw_sens)
        for i=1:numel(imus)
            yaw_mes = SensorMeasurement(min_samples,1);
            [yaw_mes.type] = deal(yaw_sens.binary_type);
            [yaw_mes.size] = deal(numel(yaw_sens.measurement));
            data = zeros(min_samples,1);
            yaw_mes.setMesArray(data);
            mes = [mes yaw_mes];
            
            labels{end+1} = [imus(i).name '_YAW'];
        end
    end

    %Note, this function assumes IMU data is already aligned and
    %ready to go. It will use the imu with least samples for the
    %output
end