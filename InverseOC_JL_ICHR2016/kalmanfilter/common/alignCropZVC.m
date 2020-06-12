function [timeArray, indArray, featureSet_imu_aligned, featureSet_mocap_aligned, dataInstance_imu_aligned, dataInstance_mocap_aligned] = ...
    alignCropZVC(featureSet_mocap, featureSet_imu, dataInstance_mocap, dataInstance_imu, alignmentSignalName, allJointStrMocap, allJointStrImu, fieldsToAlign, mocap_first, mocap_last)
    % find the the zvc after the alignment point
    currJointInd_mocap = find(ismember(allJointStrMocap, alignmentSignalName), 1);
    
    targetDq = featureSet_mocap.dq(:, currJointInd_mocap);
    targetDq = filter_dualpassBW(targetDq);
    inds = calcZVC(targetDq);
    firstZVC = find(inds > mocap_first, 2, 'first');
    lastZVC = find(inds < mocap_last, 2, 'last');
    
    threshold = 25;
    
    % if the detected zvc is too close to the peak, it's probably a jitter
    if abs(mocap_first - inds(firstZVC(1))) < threshold
        firstZVC = firstZVC(2);
    else
        firstZVC = firstZVC(1);
    end
    
    % if the detected zvc is too close to the peak, it's probably a jitter
    if abs(mocap_last - inds(lastZVC(2))) < threshold
        lastZVC = lastZVC(1);
    else
        lastZVC = lastZVC(2);
    end
        
    indArray = inds(firstZVC):inds(lastZVC);
    timeArray = featureSet_mocap.time(indArray);
    
    featureSet_imu_aligned = generateCroppedDataFeatureSet(featureSet_imu, fieldsToAlign, indArray, timeArray);
    featureSet_mocap_aligned = generateCroppedDataFeatureSet(featureSet_mocap, fieldsToAlign, indArray, timeArray);
    
    dataInstance_imu_aligned = generateCroppedDataDataInstance(dataInstance_imu, fieldsToAlign, indArray, timeArray);
    dataInstance_mocap_aligned = generateCroppedDataDataInstance(dataInstance_mocap, fieldsToAlign, indArray, timeArray);
end

function inds = calcZVC(dq)
    inds = [];
    minTol = 1e-3;
    for i = 1:length(dq)-1
        win = dq(i:i+1);
        
%         if win(1) > minTol && win(2) < minTol
%             inds = [inds; i];
        if win(1) < -minTol && win(2) > -minTol
            inds = [inds; i];
        end
    end
end

function featureSet_aligned = generateCroppedDataFeatureSet(featureSet_orig, fieldToAlignNames, indArray, timeArray)
    featureSet_aligned = rlFeatureSet_ioc();
    featureSet_aligned.time = timeArray;
    featureSet_aligned.dt = featureSet_orig.dt;
    featureSet_aligned.joint_labels = featureSet_orig.joint_labels;
    featureSet_aligned.measurement_labels = featureSet_orig.measurement_labels;

    % align fields in features
    for i = 1:length(fieldToAlignNames)
        fieldToAlignStr = fieldToAlignNames{i};
        fieldToAlign = featureSet_orig.(fieldToAlignStr);
        
        switch class(fieldToAlign)
            case 'double'
                fieldToAlignDouble = fieldToAlign;
                
            case 'SensorMeasurement'
                fieldToAlignDouble = fieldToAlign.getMesArray;
        end
        
        resampledData = fieldToAlignDouble(indArray, :);
        
        switch class(fieldToAlign)
            case 'double'
                alignedData = resampledData;
                
            case 'SensorMeasurement'
                numTimeLength = size(resampledData, 1);
                numSensors = size(fieldToAlign, 2);
                alignedData = SensorMeasurement(numTimeLength, numSensors);
                
                for j = 1:numSensors
                    sensorArray = SensorMeasurement(numTimeLength, 1);
                    [sensorArray.type] = deal(fieldToAlign(1, j).type);
                    [sensorArray.size] = deal(fieldToAlign(1, j).size);
                    alignedData(:, j) = sensorArray;
                end
                
                alignedData.setMesArray(resampledData);
        end
        
        featureSet_aligned.(fieldToAlignStr) = alignedData;
    end
    
    if ~isempty(featureSet_orig.frameData)
        for i = 1:length(featureSet_orig.frameData)
            featureSet_aligned.frameData(i).position = featureSet_orig.frameData(i).position(indArray, :);
            featureSet_aligned.frameData(i).velocity = featureSet_orig.frameData(i).velocity(indArray, :);
            featureSet_aligned.frameData(i).acceleration = featureSet_orig.frameData(i).acceleration(indArray, :);
            featureSet_aligned.frameData(i).name = featureSet_orig.frameData(i).name;
        end
    end
end

function dataInstance_aligned =  generateCroppedDataDataInstance(dataInstance_orig, fieldToAlignNames, indArray, timeArray)
    if isempty(dataInstance_orig)
        dataInstance_aligned = [];
        return;
    end
    
    dataInstance_aligned = rlDataInstance(dataInstance_orig.modelInstance);
	dataInstance_aligned.dt = dataInstance_orig.dt;
    dataInstance_aligned.time = timeArray;
    
    padWithData = 1;
    
    sizeData = size(dataInstance_orig.data);
    
    if sizeData == 1
        type = 'mocap';
    else
        type = 'imu';
    end
    
    switch type
        case 'mocap'
            % probably mocap
            fieldNames = fieldnames(dataInstance_orig.data);
            for i = 1:length(fieldNames)
                fieldToAlignStr = fieldNames{i};
                fieldToAlign = dataInstance_orig.data.(fieldToAlignStr);
                alignedData = fieldToAlign(indArray, :);
                dataInstance_aligned.data.(fieldToAlignStr) = alignedData;
            end
            
        case 'imu'
            % probably imu
            for i = 1:length(dataInstance_orig.data)
                fieldNames = fieldnames(dataInstance_orig.data(i));
                
                for j = 1:length(fieldNames)
                    fieldToAlignStr = fieldNames{j};
                    fieldToAlign = dataInstance_orig.data(i).(fieldToAlignStr);
                    
                    switch fieldToAlignStr
                        case {'accelerometerCalibrated', 'gyroscopeCalibrated'}
                            alignedData = fieldToAlign(indArray, :);
                            
                        otherwise
                            alignedData = fieldToAlign;
                    end
                    
                    dataInstance_aligned.data(i).(fieldToAlignStr) = alignedData;
                end
            end
    end
     
    dataInstance_aligned.measurement_labels = dataInstance_orig.measurement_labels;
	dataInstance_aligned.measurement_sensor = [];
	dataInstance_aligned.segments = [];
end