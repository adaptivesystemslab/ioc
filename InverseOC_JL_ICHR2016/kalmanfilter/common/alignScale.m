function [dt, timeArray, featureSet_imu_aligned, featureSet_mocap_aligned, dataInstance_mocap_aligned, dataInstance_imu_aligned, h, ...
    imu_first, imu_last, mocap_first, mocap_last, currJointInd_imu, currJointInd_mocap, ...
    startingOffset, endingOffset, output_N] = alignScale(featureSet_mocap, featureSet_imu, dataInstance_mocap, dataInstance_imu, ...
    alignmentSignalName, allJointStrMocap, allJointStrImu, fieldsToAlign)
    % THE SETTINGS WE MAY HAVE TO CHANGE ARE HERE
    %These are the threshold for the IMUs
    startingOffset = 100;
    endingOffset = 50;
%     thresholds_first = 0.002;
%     thresholds_last = 0.002; % s_square
    thresholds_first = 1;
    thresholds_last = 1; % s_square

    currJointInd_mocap = find(ismember(allJointStrMocap, alignmentSignalName), 1);
    currJointInd_imu = find(ismember(allJointStrImu, alignmentSignalName), 1);
    
    [imu_first, imu_last, imu_alignmentSignal] = alignment_q_startEnd(featureSet_imu, currJointInd_imu, thresholds_first, thresholds_last);
    [mocap_first, mocap_last, mocap_alignmentSignal] = alignment_q_startEnd(featureSet_mocap, currJointInd_mocap, thresholds_first, thresholds_last);
    
    %Get the mean number of points between imu peaks
%     output_N = floor(mean(imu_last-imu_first));
% startingOffset = imu_first - 1;
% endingOffset = size(featureSet_imu.q, 1) - imu_last;
    output_N = floor(mean(mocap_last-mocap_first));
    startingOffset = mocap_first - 1;
    endingOffset = size(featureSet_mocap.q, 1) - mocap_last;

    imu_delta = floor(mean(imu_last-imu_first));
    mocap_delta = floor(mean(mocap_last-mocap_first));
    imu_indexes = alignment_index(imu_first, imu_last, output_N, startingOffset, endingOffset);
    mocap_indexes = alignment_index(mocap_first, mocap_last, output_N, startingOffset, endingOffset);

    dt = mean(diff(featureSet_mocap.time));
    timeArray = (0:length(mocap_indexes)-1)*dt;
    
    featureSet_imu_aligned = generateAlignedDataFeatureSet(featureSet_imu, output_N, imu_delta, imu_indexes, fieldsToAlign, timeArray, dt);
    featureSet_mocap_aligned = generateAlignedDataFeatureSet(featureSet_mocap, output_N, mocap_delta, mocap_indexes, fieldsToAlign, timeArray, dt);
    
    dataInstance_imu_aligned = generateAlignedDataDataInstance(dataInstance_imu, output_N, imu_delta, imu_indexes, fieldsToAlign, timeArray, dt);
    dataInstance_mocap_aligned = generateAlignedDataDataInstance(dataInstance_mocap, output_N, mocap_delta, mocap_indexes, fieldsToAlign, timeArray, dt);
    
%     q_imu = alignment_resampleAndIndex(featureSet_imu.q, imu_N, imu_delta, imu_indexes);
%     dq_imu = alignment_resampleAndIndex(featureSet_imu.dq, imu_N, imu_delta, imu_indexes);
%     x_imu = alignment_resampleAndIndex(featureSet_imu.x, imu_N, imu_delta, imu_indexes);
%     q_mocap = alignment_resampleAndIndex(featureSet_mocap.q, imu_N, mocap_delta, mocap_indexes);
% %     dq_mocap = calcDerivVert(q_mocap, featureSet_imu.dt);
%     dq_mocap = alignment_resampleAndIndex(featureSet_mocap.dq, imu_N, mocap_delta, mocap_indexes);
%     x_mocap = alignment_resampleAndIndex(featureSet_mocap.x, imu_N, mocap_delta, mocap_indexes);
    
%     if 0
%         h = figure; 
%         time_orig = 1:length(featureSet_imu.q);
%         q_orig = featureSet_imu.q;
%         time_align = imu_indexes;
%         q_align = q_imu;
%         
%         plot(time_orig, q_orig);
%         hold on;
%         plot(time_align, q_align, '--');
%     end
    
    if 0
        h = figure;
        signal_imu = filter_dualpassBW(featureSet_imu.q);
        signal_mocap = filter_dualpassBW(featureSet_mocap.q);
        subplot(211); hold on;
        plot(signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
        plot(imu_first, signal_imu(imu_first, currJointInd_imu), 'kx');
        plot(imu_last, signal_imu(imu_last, currJointInd_imu), 'kx');
        plot(signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap', 'LineStyle', '--');
        plot(mocap_first, signal_mocap(mocap_first, currJointInd_mocap), 'kx');
        plot(mocap_last, signal_mocap(mocap_last, currJointInd_mocap), 'kx');
        title('orig filt(q)');
        
        %         figure;
        signal_imu = featureSet_imu_aligned.q;
        signal_mocap = featureSet_mocap_aligned.q;
        subplot(212);  hold on;
        plot(signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
        plot(startingOffset, signal_imu(startingOffset, currJointInd_imu), 'kx');
        plot(size(signal_imu, 1)-endingOffset, signal_imu(size(signal_imu, 1)-endingOffset, currJointInd_imu), 'kx');
        plot(signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap', 'LineStyle', '--');
        plot(startingOffset, signal_mocap(startingOffset, currJointInd_mocap), 'kx');
        plot(size(signal_imu, 1)-endingOffset, signal_mocap(size(signal_imu, 1)-endingOffset, currJointInd_mocap), 'kx');
        title('output q');
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

function [firstInd, lastInd, alignmentSignal] = alignment_q_startEnd(featureSet, currJointInd, thresholds_first, thresholds_last)
    alignmentSignal = featureSet.q(:, currJointInd);
    alignmentSignal = normVector(alignmentSignal);
%     alignmentSignal = filter_dualpassBW(alignmentSignal);
    
    % extract the peak from s_square
    [imu_peaks,imu_locs] = findpeaks(alignmentSignal);
    firstInd = imu_locs(find(alignmentSignal(imu_locs)>thresholds_first,1));
    lastInd = imu_locs(find(alignmentSignal(imu_locs)>thresholds_last,1,'last'));
end

function [firstInd, lastInd] = alignment_dq_startEnd(featureSet, currJointInd, half_w_d, thresholds_first, thresholds_last)
    alignmentSignal = featureSet.dq(:, currJointInd);
    alignmentSignal = filter_dualpassBW(alignmentSignal);
    
    total_sample = size(alignmentSignal,1);
    normal_x = zeros(1, total_sample);
    for i = (half_w_d + 1):total_sample - (half_w_d + 1)
        normal_x(i) = var(alignmentSignal(i - half_w_d:i + half_w_d,1));
    end
    s_square = normal_x.^2;
    
    % extract the peak from s_square
    [imu_peaks,imu_locs] = findpeaks(s_square);
    imu_first = imu_locs(find(s_square(imu_locs)>thresholds_first,1));
    imu_last = imu_locs(find(s_square(imu_locs)>thresholds_last,1,'last'));
    
    % now find the peak from the original
    normSignal = normVector(alignmentSignal);
    [norm_peaks, norm_locs] = findpeaks(normSignal);
    norm_first = norm_locs(find(norm_locs>imu_first,2));
    norm_last = norm_locs(find(norm_locs<imu_last,2,'last'));
    
    % figure out which points to align on. 
    % starting dq should be neg than pos
    if alignmentSignal(norm_first(1)) < 0 && alignmentSignal(norm_first(2)) > 0
        firstInd = norm_first;
    else
        % first peak is missing
        firstInd(1) = 0;
        firstInd(2) = norm_first(1);
    end
    
    % ending dq should be pos than neg
    if alignmentSignal(norm_last(1)) < 0 && alignmentSignal(norm_last(2)) > 0
        lastInd = norm_last;
    else
        % last peak is missing
        lastInd(1) = norm_last(2);
        lastInd(2) = 0;
    end
    
    % lastly, go back to the original position plot
    alignmentSignal = featureSet.q(:, currJointInd);
    [minVal, firstInd] = max(abs(alignmentSignal(1:firstInd(2))));
    [minVal, lastInd_tmp] = max(abs(alignmentSignal(lastInd(1):end)));
    lastInd = lastInd_tmp + lastInd(1);
    
%     firstInd = norm_first;
%     lastInd = norm_last;
end

function indexes = alignment_index(first, last, N, startingOffset, endingOffset)
    first_out = first*N/(last-first);
    last_out = last*N/(last-first);
    indexes = round(first_out-startingOffset:last_out+endingOffset);
end

function featureSet_aligned = generateAlignedDataFeatureSet(featureSet_orig, imu_N, imu_delta, imu_indexes, fieldToAlignNames, timeArray, dt)
    featureSet_aligned = rlFeatureSet_ioc();
    featureSet_aligned.time = timeArray;
    featureSet_aligned.dt = dt;
    featureSet_aligned.joint_labels = featureSet_orig.joint_labels;
    featureSet_aligned.measurement_labels = featureSet_orig.measurement_labels;
    
    padWithData = 1;
    
    % align fields in features
    for i = 1:length(fieldToAlignNames)
        fieldToAlignStr = fieldToAlignNames{i};
        fieldToAlign = featureSet_orig.(fieldToAlignStr);
        
        if isempty(fieldToAlign)
            featureSet_aligned.(fieldToAlignStr) = [];
            continue;
        end

%         switch fieldsToAlign{i}
% %             case {'dq', 'baseVelocity', 'baseAcceleration'}
% %                 padWithData = 0;
%                 
%             otherwise
%                 padWithData = 1;
%         end

        switch class(fieldToAlign)
            case 'double'
                alignedData = alignment_resampleAndIndex(fieldToAlign, imu_N, imu_delta, imu_indexes, padWithData);
        
            case 'SensorMeasurement'
                fieldToAlignDouble = fieldToAlign.getMesArray;
                resampledData = alignment_resampleAndIndex(fieldToAlignDouble, imu_N, imu_delta, imu_indexes, padWithData);
                
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
    
    % also align fields in frameData
    if ~isempty(featureSet_orig.frameData)
        for i = 1:length(featureSet_orig.frameData)
            featureSet_aligned.frameData(i).position = alignment_resampleAndIndex(featureSet_orig.frameData(i).position, imu_N, imu_delta, imu_indexes, padWithData);
            featureSet_aligned.frameData(i).velocity = alignment_resampleAndIndex(featureSet_orig.frameData(i).velocity, imu_N, imu_delta, imu_indexes, padWithData);
            featureSet_aligned.frameData(i).acceleration = alignment_resampleAndIndex(featureSet_orig.frameData(i).acceleration, imu_N, imu_delta, imu_indexes, padWithData);
            featureSet_aligned.frameData(i).name = featureSet_orig.frameData(i).name;
        end
    end
end

function dataInstance_aligned = generateAlignedDataDataInstance(dataInstance_orig, imu_N, imu_delta, imu_indexes, fieldToAlignNames, timeArray, dt)
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
                alignedData = alignment_resampleAndIndex(fieldToAlign, imu_N, imu_delta, imu_indexes, padWithData);
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
                            alignedData = alignment_resampleAndIndex(fieldToAlign, imu_N, imu_delta, imu_indexes, padWithData);
                            
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

