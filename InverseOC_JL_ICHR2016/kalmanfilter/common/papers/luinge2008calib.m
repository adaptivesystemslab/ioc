function transform = luinge2008calib(dataInstance, modelInstance)
    % https://ieeexplore.ieee.org/document/1273529
    % determine the direction of gravity in all 5 sensors
    gravInds = 1:10;
    
    for i = 1:length(dataInstance.data)
        currImuData = dataInstance.data(i);
        
        % use gravity to get Z
        accelData = currImuData.accelerometerCalibrated(gravInds, :);
        meanAccelData = mean(accelData, 1);
        currZVec = -meanAccelData / norm(meanAccelData);
        
        % use hip ext to get Y
        gyroData = currImuData.gyroscopeCalibrated;
        normGyroData = abs(normVector(currImuData.gyroscopeCalibrated));
        [maxVal, maxInd] = max(normGyroData);
        currYVec = gyroData(maxInd, :) / norm(gyroData(maxInd, :));
        
        currXVec = cross(currYVec, currZVec);
        currZVec = cross(currXVec, currYVec);
        
        currRotMtx = [currXVec' currYVec' currZVec'];
        
        % proper normalizing
        [S, V, D] = svd(currRotMtx);
        V = eye(3);
        currRotMtx = S*V*D';

        sensorStr = currImuData.name;
        frame.name = modelInstance.sensorSecondaryAttachmentArray{i}{1};
        
        switch currImuData.name
            case 'HIP_BASE_JC'
                T_f_imu = eye(4);
                T_f_P = eye(4);
                T_f_Q = eye(4);
                T_f_R = eye(4);
                
            otherwise
                T_f_imu = currRotMtx;
                T_f_P = eye(4);
                T_f_Q = eye(4);
                T_f_R = eye(4);
        end
        
        transform(i) = imuTrans(sensorStr{2}, frame, T_f_imu, T_f_P, T_f_Q, T_f_R);
    end
end