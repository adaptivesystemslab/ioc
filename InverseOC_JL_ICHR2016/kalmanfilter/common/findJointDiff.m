function dist = findJointDiff(featureSet_mocap, featureSet_imu, frameToDiff_imu, frameToDiff_mocap)
    % find the mocap frames
    allMocapNames = featureSet_mocap.joint_labels;
    allImuNames = featureSet_imu.joint_labels;
    
    q_mocap = rad2deg(featureSet_mocap.q);
    q_imu = rad2deg(featureSet_imu.q);
    
    rmseFct = setRmseFct();
    
    if iscell(frameToDiff_mocap)
        for i = 1:length(frameToDiff_mocap)
            indMocapFrame = find(ismember(allMocapNames, frameToDiff_mocap{i}) == 1);
            indImuFrame = find(ismember(allImuNames, frameToDiff_imu{i}) == 1);
            
            currRmse = rmseFct(q_mocap(:, indMocapFrame), q_imu(:, indImuFrame));
            dist(i) = currRmse;
        end
    else
        indMocapFrame = find(ismember(allMocapNames, frameToDiff_mocap) == 1);
        indImuFrame = find(ismember(allImuNames, frameToDiff_imu) == 1);
        
        currRmse = rmseFct(q_mocap(:, indMocapFrame), q_imu(:, indImuFrame));
        dist = currRmse;
    end
end