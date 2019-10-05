function trc = loadDataFromTrc_healthy1_arsLoader(pathToRawHeader, pathToRawData, chainDirection)
    mocapData = mocapDataHandle_crop(pathToRawHeader, pathToRawData);
    mocapData.load();
    
    dt = 0.01;
    startTime = mocapData.time(1);
    endTime = mocapData.time(end);
    time = startTime:dt:endTime;
 
    if length(time) < 5
        error('Time length malformed');
    end
    
    field_names = fieldnames(mocapData.data);
    field_names = field_names(4:end);
    
    % interpolate data to be evenly arrayed
    for i=1:numel(field_names)
        mocapData.data.(field_names{i}) = interp1(mocapData.time, mocapData.data.(field_names{i}), time)' / 1000; % convert to meters
    end
    
    trc.dt = dt;
    trc.time = time;
    
    % construct XYZ of markers
    trc.data.SHOULDER_R  = [mocapData.data.R_Shoulder_x      mocapData.data.R_Shoulder_y      mocapData.data.R_Shoulder_z];
    trc.data.SHOULDER_L  = [mocapData.data.L_Shoulder_x      mocapData.data.L_Shoulder_y      mocapData.data.L_Shoulder_z];
    trc.data.ASIS_R      = [mocapData.data.R_ASIS_x          mocapData.data.R_ASIS_y          mocapData.data.R_ASIS_z];
    trc.data.ASIS_L      = [mocapData.data.L_ASIS_x          mocapData.data.L_ASIS_y          mocapData.data.L_ASIS_z];
    trc.data.KNEE_R_MED  = [mocapData.data.R_Knee_Medial_x   mocapData.data.R_Knee_Medial_y   mocapData.data.R_Knee_Medial_z];
    trc.data.KNEE_R_LAT  = [mocapData.data.R_Knee_Lateral_x  mocapData.data.R_Knee_Lateral_y  mocapData.data.R_Knee_Lateral_z];
    trc.data.ANKLE_R_MED = [mocapData.data.R_Ankle_Medial_x  mocapData.data.R_Ankle_Medial_y  mocapData.data.R_Ankle_Medial_z];
    trc.data.ANKLE_R_LAT = [mocapData.data.R_Ankle_Lateral_x mocapData.data.R_Ankle_Lateral_y mocapData.data.R_Ankle_Lateral_z];
%     trc.data.TOE_R       = [mocapData.data.R_Toe_x           mocapData.data.R_Toe_y           mocapData.data.R_Toe_z];
%     trc.data.HEEL_R      = [mocapData.data.R_Heel_x          mocapData.data.R_Heel_y          mocapData.data.R_Heel_z];
    
    trc.data.HIP_R_IMU   = [mocapData.data.Hip_Shimmer_x     mocapData.data.Hip_Shimmer_y     mocapData.data.Hip_Shimmer_z];
    trc.data.KNEE_R_IMU  = [mocapData.data.Knee_Shimmer_x    mocapData.data.Knee_Shimmer_y    mocapData.data.Knee_Shimmer_z];
    trc.data.ANKLE_R_IMU = [mocapData.data.Ankle_Shimmer_x   mocapData.data.Ankle_Shimmer_y   mocapData.data.Ankle_Shimmer_z];
    
    if 0
        vis = rlVisualizer('vis',640,480);

        field_names = fieldnames(trc.data);
        for j = 1:length(trc.time)
            for i=1:numel(field_names)
                vis.addMarker(field_names{i}, trc.data.(field_names{i})(j, 1:3));
            end
            pause(0.01);
            vis.update();
        end
    end

%     % construct averaged position of markers TODO
%     trc.data.TOE_JC_R = [trc.data.TOE_R];
%     trc.data.ANKLE_JC_R = [trc.data.ANKLE_R_LAT + trc.data.ANKLE_R_MED] / 2;
%     trc.data.KNEE_JC_R = [trc.data.KNEE_R_LAT + trc.data.KNEE_R_MED] / 2;
%     trc.data.HIP_JC_R = [trc.data.ASIS_R];
%     trc.data.SHOULDER_JC_R = [trc.data.SHOULDER_R];
    
% %     % remove offset
%     switch chainDirection
%         case 'ankle'
%             offsetVal = trc.data.ANKLE_JC_R;
%         
%         case 'hip'
%             offsetVal = trc.data.HIP_JC_R;
%             
%         case 'shoulder'
%             offsetVal = trc.data.SHOULDER_JC_R;
%     end
    
    field_names = fieldnames(trc.data);

    for i=1:numel(field_names)
%         trc.data.(field_names{i}) = trc.data.(field_names{i}) - offsetVal;
        
%         rotate data so it's in X and Y
%          trc.data.(field_names{i}) = (rotz(pi/2)*trc.data.(field_names{i})')';
%         trc.data.(field_names{i}) = (rotz(ang1)*trc.data.(field_names{i})')';
        
%         % project to 2D
%         trc.data.(field_names{i}) = [ ...
%                 zeros(size(trc.data.(field_names{i})(:, 1))) ...
%                 trc.data.(field_names{i})(:, 2) ...
%                 trc.data.(field_names{i})(:, 3) ...
%                 ];
    end
    
    if 0
        vis = rlVisualizer('vis',640,480);

        field_names = fieldnames(trc.data);
        for j = 1
            for i=1:numel(field_names)
                vis.addMarker(field_names{i}, trc.data.(field_names{i})(j, 1:3));
            end
            pause(0.01);
            vis.update();
        end
    end
end