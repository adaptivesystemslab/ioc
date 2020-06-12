function [dt, timeArray, featureSet_imu_aligned, featureSet_mocap_aligned, dataInstance_mocap_aligned, dataInstance_imu_aligned, h, ...
    D] = alignXcorr(featureSet_mocap, featureSet_imu, dataInstance_mocap, dataInstance_imu, alignmentSignalName, allJointStrMocap, allJointStrImu, fieldsToAlignNames)
    % THE SETTINGS WE MAY HAVE TO CHANGE ARE HERE
    %These are the threshold for the IMUs
    maxlag = 100;
    
    currJointInd_mocap = find(ismember(allJointStrMocap, alignmentSignalName), 1);
    currJointInd_imu = find(ismember(allJointStrImu, alignmentSignalName), 1);
        
    signal_mocap = featureSet_mocap.q(:, currJointInd_mocap);
    signal_imu = featureSet_imu.q(:, currJointInd_imu);
    
    [signal_mocap_align, signal_imu_align, D] = alignsignals(signal_mocap,signal_imu,maxlag,'truncate');

    featureSet_mocap_aligned = featureSet_mocap; % don't need to change the mocap
    
    featureSet_imu_aligned = generateAlignedDataFeatureSet(signal_imu, fieldsToAlignNames, featureSet_imu, D);
    
    dt = featureSet_mocap_aligned.dt;
    timeArray = featureSet_mocap_aligned.time;
    
    featureSet_imu_aligned.dt = featureSet_mocap_aligned.dt;
    featureSet_imu_aligned.time = featureSet_mocap_aligned.time;
    featureSet_imu_aligned.joint_labels = featureSet_imu.joint_labels;
    featureSet_imu_aligned.measurement_labels = featureSet_imu.measurement_labels;
    
    lenSig = size(signal_imu, 1);
    dataInstance_mocap_aligned = dataInstance_mocap; % don't need to change the mocap
    
    if isempty(dataInstance_imu)
        dataInstance_imu_aligned = [];
        h = [];
        return;
    end
    
    dataInstance_imu_aligned = rlDataInstance_imu(dataInstance_imu.modelInstance);
    
    dataInstance_imu_aligned.dt = featureSet_mocap_aligned.dt;
    dataInstance_imu_aligned.time = timeArray;
    
    for i = 1:length(dataInstance_imu.data)
        fieldNames = fieldnames(dataInstance_imu.data(i));
        
        for j = 1:length(fieldNames)
            fieldToAlignStr = fieldNames{j};
            fieldToAlignDouble = dataInstance_imu.data(i).(fieldToAlignStr);
            
            switch fieldToAlignStr
                case {'accelerometerCalibrated', 'gyroscopeCalibrated'}
                    if D > 0
                        padVals = repmat(fieldToAlignDouble(lenSig, :), D, 1);
                        resampledData = [ fieldToAlignDouble((D+1):lenSig, :); ...
                            padVals;];
                    elseif D < 0
                        padVals = repmat(fieldToAlignDouble(1, :), abs(D), 1);
                        resampledData = [padVals; ...
                            fieldToAlignDouble(1:lenSig-abs(D), :)];
                    else
                        resampledData = fieldToAlignDouble;
                    end
                    
                    alignedData = resampledData;
                    dataInstance_imu_aligned.data(i).(fieldToAlignStr) = alignedData;
            end
        end
    end
    
    dataInstance_imu_aligned.measurement_labels = dataInstance_imu.measurement_labels;
    dataInstance_imu_aligned.measurement_sensor = [];
    dataInstance_imu_aligned.segments = [];
    
    if 0
        h = figure;
        signal_imu = filter_dualpassBW(featureSet_imu.q);
        signal_mocap = filter_dualpassBW(featureSet_mocap.q);
        subplot(211); hold on;
        plot(signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
        
        plot(signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap');
        title('orig (q)');
        
%         figure;
        signal_imu = featureSet_imu_aligned.q;
        signal_mocap = featureSet_mocap_aligned.q;
        subplot(212);  hold on;
        plot(signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
        plot(signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap');
        title(['output q (D = ' num2str(D) ')']);
    else 
        h = [];
    end
end

function signal_out = alignment_resampleAndIndex(signal_orig, imu_N, delta, indexes, prepadWithData)
    % if the indexes are out of bounds, pre/postpend entries
    if min(indexes) < 1 
        shift = abs(min(indexes)) + 1;
        
        signal_toAugment = signal_orig(1, :);
        if ~prepadWithData
            signal_toAugment = zeros(size(signal_toAugment));            
        end
        signal_augment = repmat(signal_toAugment, shift, 1);
        
        signal_orig = [signal_augment; signal_orig];
        indexes = indexes + shift;
    end
    
    if max(indexes) > size(signal_orig, 2)
        shift = max(indexes) - size(signal_orig, 2);
        
        signal_toAugment = signal_orig(end, :);
        if ~prepadWithData
            signal_toAugment = zeros(size(signal_toAugment));
        end
        signal_augment = repmat(signal_toAugment, shift, 1);
        
        signal_orig = [signal_orig; signal_augment];
    end

    signal_resampled = resample(signal_orig,imu_N,delta,0);
    signal_out = signal_resampled(indexes, :);
end

function featureSet_imu_aligned = generateAlignedDataFeatureSet(signal_imu, fieldsToAlignNames, featureSet_imu, D)
    featureSet_imu_aligned = rlFeatureSet_ioc();
    
    lenSig = size(signal_imu, 1);
    for i = 1:length(fieldsToAlignNames)
        fieldToAlignStr = fieldsToAlignNames{i};
        fieldToAlign = featureSet_imu.(fieldToAlignStr);
        
        switch class(fieldToAlign)
            case 'double'
                fieldToAlignDouble = fieldToAlign;
                
            case 'SensorMeasurement'
                fieldToAlignDouble = fieldToAlign.getMesArray;
        end
        
        if D > 0
            padVals = repmat(fieldToAlignDouble(lenSig, :), D, 1);
            resampledData = [ fieldToAlignDouble((D+1):lenSig, :); ...
                padVals;];
        elseif D < 0
            padVals = repmat(fieldToAlignDouble(1, :), abs(D), 1);
            resampledData = [padVals; ...
                fieldToAlignDouble(1:lenSig-abs(D), :)];
        else
            resampledData = fieldToAlignDouble;
        end
        
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
        
        featureSet_imu_aligned.(fieldToAlignStr) = alignedData;
    end
end