function rmseVal = findFrameCartRmse(featureSet_mocap, featureSet_imu, frameToDiff)
    % find the mocap frames
    allMocapNames = {featureSet_mocap.frameData.name};
    allImuNames = {featureSet_imu.frameData.name};
    indMocapFrame = find(ismember(allMocapNames, frameToDiff{1}) == 1);
    indImuFrame = find(ismember(allImuNames, frameToDiff{1}) == 1);
    
    % calculate the rotation between frame 1 and 2
    dist = zeros(length(featureSet_mocap.time), 1);
    for i = 1:length(featureSet_mocap.time)
        R_mocap1 = reshape(featureSet_mocap.frameData(indMocapFrame).position(i, :), 4, 4);
        mocapPos(i, :) = R_mocap1(1:3, 4);
        
        R_imu1 = reshape(featureSet_imu.frameData(indImuFrame).position(i, :), 4, 4);
        imuPos(i, :) = R_imu1(1:3, 4);
    end
    
    rmseFct = setRmseFct();
    rmseVal = rmseFct(imuPos, mocapPos);
end