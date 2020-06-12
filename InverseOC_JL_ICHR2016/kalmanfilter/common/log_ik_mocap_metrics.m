function [markerRmse, jointRom, jointMaxDq, jointMaxDdq, h3] = calc_ik_mocap_metrics(modelInstance, featureSet_mocap,...
    ekf_markerMask, ekf_eventhandler, ekf_markerTally, filepathMocapIkLog ,currFileEntry, externalParam)
    
    % calculate mocap IK metrics  
    [jointRom, jointMaxDq, jointMaxDdq] = calc_rom_dq_ddq(featureSet_mocap);
    
%     [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.swappedmissing);
%     markerUsedPercent = sum(markerUsed)/sum(markerTotal);
    
    %         write header to file
    if ~exist(filepathMocapIkLog, 'file')
        %             if doesn't exist write header
        header = 'Subject,Exercise,linkDef,suffix';
        header = [header ',procNoiseMocapPris,procNoiseMocapRev,markerSwapThres'];
        header = [header ',obsNoiseMarkerPos,obsNoiseMarkerVel'];
        header = [header ',marker_total,marker_correct,marker_swapped,marker_missing'];
        header = [header ',marker_correct_percent,marker_swapped_percent,marker_missing_percent'];
        
        header = [header ',maskSwapMiss_markerUsed_sum,maskSwapMiss_markerUsed_mean,maskSwapMiss_markerUsed_stddev,maskSwapMiss_markerRMSE_mean [m],maskSwapMiss_markerRMSE_stddev'];
        header = [header ',maskNone_markerUsed_sum,maskNone_markerUsed_mean,maskNone_markerUsed_stddev,maskNone_markerRMSE_mean [m],maskNone_markerRMSE_stddev'];
        header = [header ',rom_mean,rom_stddev'];
        header = [header ',maxdq_mean,maxdq_stddev'];
        
        prefix = 'maskSwapMiss';
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerUsed_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerRMSE_mean_' featureSet_mocap.measurement_labels{j}];
        end
        
        prefix = 'maskNone';
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerUsed_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerRMSE_mean_' featureSet_mocap.measurement_labels{j}];
        end
        
        prefix = 'maskSwap';
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerUsed_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerRMSE_mean_' featureSet_mocap.measurement_labels{j}];
        end
        
        prefix = 'maskMiss';
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerUsed_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' prefix '_markerRMSE_mean_' featureSet_mocap.measurement_labels{j}];
        end
        
        for i = 1:length(featureSet_mocap.joint_labels)
            header = [header ',rom_' featureSet_mocap.joint_labels{i}];
        end
        for i = 1:length(featureSet_mocap.joint_labels)
            header = [header ',maxDq_' featureSet_mocap.joint_labels{i}];
        end
        %         for i = 1:length(featureSet_mocap.joint_labels)
        %             header = [header ',maxDdq_' featureSet_mocap.joint_labels{i}];
        %         end
        
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' 'marker_total_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' 'marker_correct_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' 'marker_swapped_' featureSet_mocap.measurement_labels{j}];
        end
        for j = 1:length(featureSet_mocap.measurement_labels)
            header = [header ',' 'marker_missing_' featureSet_mocap.measurement_labels{j}];
        end
        
        header = [header '\n'];
    else
        header = '';
    end

    % open/create log file
    [fileID, errmsg] = fopen(filepathMocapIkLog, 'a');
    fprintf(fileID, header);
    
    % write data to log file
    fprintf(fileID,'%s,%s,%s,%s',currFileEntry.subjectString, currFileEntry.exerciseName, ...
        externalParam.linkDefinition, externalParam.subSuffix);
    
    if isfield(featureSet_mocap.ekfParams, 'processNoiseCoefficientMocapPrismatic')
        fprintf(fileID,',%f,%f,%f',featureSet_mocap.ekfParams.processNoiseCoefficientMocapPrismatic, ...
            featureSet_mocap.ekfParams.processNoiseCoefficientMocapRevolute, ...
            featureSet_mocap.ekfParams.markerSwapTolerance);
    else
        fprintf(fileID,',%f,%f,%f',featureSet_mocap.ekfParams.processNoiseCoefficient, ...
            0, ...
            featureSet_mocap.ekfParams.markerSwapTolerance);
    end
    
    fprintf(fileID,',%f,%f', ...
        featureSet_mocap.ekfParams.observationNoiseCoefficientMarkerPosition, ...
        featureSet_mocap.ekfParams.observationNoiseCoefficientMarkerVelocity);
   
    fprintf(fileID,',%u,%u,%u,%u', sum(ekf_markerTally.total), ...
        sum(ekf_markerTally.correct), ...
        sum(ekf_markerTally.swapped), ...
        sum(ekf_markerTally.missing));
    fprintf(fileID,',%f,%f,%f', ...
        sum(ekf_markerTally.correct)/sum(ekf_markerTally.total), ...
        sum(ekf_markerTally.swapped)/sum(ekf_markerTally.total), ...
        sum(ekf_markerTally.missing)/sum(ekf_markerTally.total));
    
    [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.swappedmissing);
    if ~isempty(featureSet_mocap.measurement_output)
        fprintf(fileID,',%.4f', sum(markerUsed)/sum(markerTotal));
        fprintf(fileID,',%.4f', mean(markerUsed./markerTotal));
        fprintf(fileID,',%.4f', std(markerUsed./markerTotal));
%         fprintf(fileID,',%.4f', mean(markerRmse(modelInstance.inds_mocapMarkers)));
%         fprintf(fileID,',%.4f', std(markerRmse(modelInstance.inds_mocapMarkers)));
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
    else
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
    end
    
    [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.none);
    if ~isempty(featureSet_mocap.measurement_output)
        fprintf(fileID,',%.4f', sum(markerUsed)/sum(markerTotal));
        fprintf(fileID,',%.4f', mean(markerUsed./markerTotal));
        fprintf(fileID,',%.4f', std(markerUsed./markerTotal));
%         fprintf(fileID,',%.4f', mean(markerRmse(modelInstance.inds_mocapMarkers)));
%         fprintf(fileID,',%.4f', std(markerRmse(modelInstance.inds_mocapMarkers)));
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
    else
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
        fprintf(fileID,',%.4f', 0);
    end
    
%     fprintf(fileID,',%.4f', mean(jointRom(modelInstance.inds_joints)));
%     fprintf(fileID,',%.4f', std(jointRom(modelInstance.inds_joints)));
%     
%     fprintf(fileID,',%.4f', mean(jointMaxDq(modelInstance.inds_joints)));
%     fprintf(fileID,',%.4f', std(jointMaxDq(modelInstance.inds_joints)));
    fprintf(fileID,',0,0,0,0');
    
    [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.swappedmissing);
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerUsed(j)/markerTotal(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerRmse(j));
    end
    
    [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.none);
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerUsed(j)/markerTotal(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerRmse(j));
    end
    
    [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.swapped);
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerUsed(j)/markerTotal(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerRmse(j));
    end
    
    [markerRmse, markerUsed, markerTotal] = calc_rmse_mocap(featureSet_mocap, ekf_markerMask.missing);
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerUsed(j)/markerTotal(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', markerRmse(j));
    end

    for j = 1:length(featureSet_mocap.joint_labels)
        fprintf(fileID,',%.4f', jointRom(j));
    end
    
    for j = 1:length(featureSet_mocap.joint_labels)
        fprintf(fileID,',%.4f', jointMaxDq(j));
    end
    
%     for j = 1:length(featureSet_mocap.joint_labels)
%         fprintf(fileID,',%.4f', jointMaxDdq(j));
%     end
    
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', ekf_markerTally.total(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', ekf_markerTally.correct(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', ekf_markerTally.swapped(j));
    end
    for j = 1:length(featureSet_mocap.measurement_labels)
        fprintf(fileID,',%.4f', ekf_markerTally.missing(j));
    end
    
    fprintf(fileID, '\n');
    
    fclose(fileID);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    figPos = [1.9282e+03 -159 1.8664e+03 932.8000];
    
    mdl = modelInstance.model;
    mdlAllJoints = {mdl.joints.name}';
    
    errorJoints_RLDefault{1} = {'frame_6dof_prismx0_to_frame_6dof_prismx1', ...
        'frame_6dof_prismy0_to_frame_6dof_prismy1', ...
        'frame_6dof_prismz0_to_frame_6dof_prismz1'};
    errorJoints_RLDefault{2} = {'joint_rhip_0', ...
        'joint_rhip_1', ...
        'joint_rhip_2'};
    errorJoints_RLDefault{3} = {'joint_rknee_0'};
    errorJoints_RLDefault{4} = {'joint_rankle_0', ...
        'joint_rankle_1'};
    
    errorJoints_RLDefault{5} = {'frame_6dof_revx0_to_frame_6dof_revx1', ...
        'frame_6dof_revy0_to_frame_6dof_revy1', ...
        'frame_6dof_revz0_to_frame_6dof_revz1'};
    errorJoints_RLDefault{6} = {'joint_lhip_0', ...
        'joint_lhip_1', ...
        'joint_lhip_2'};
    errorJoints_RLDefault{7} = {'joint_lknee_0'};
    errorJoints_RLDefault{8} = {'joint_lankle_0', ...
        'joint_lankle_1'};

    frameDofs = length(errorJoints_RLDefault);
    
    h3 = figure('position', figPos);
    ax_frame = [];
    for j = 1:frameDofs
        jointInds = [];
        for k = 1:length(errorJoints_RLDefault{j})
            jointInd = find(ismember(mdlAllJoints, errorJoints_RLDefault{j}{k}));
            if ~isempty(jointInd)
                jointInds = [jointInds; jointInd];
            end
        end
        
        mocap_time = featureSet_mocap.time;
        mocap_sig = rad2deg(featureSet_mocap.q(:, jointInds));
%         imu_sig = rad2deg(featureSet_imu.q(:, jointInds));
        
        switch length(jointInds)
            case 1
                blanks = zeros(size(mocap_sig, 1), 1);
                mocap_sig = [mocap_sig blanks blanks];
%                 imu_sig = [imu_sig blanks blanks];
                
            case 2
                blanks = zeros(size(mocap_sig, 1), 1);
                mocap_sig = [mocap_sig(:, 1) blanks mocap_sig(:, 2)];
%                 imu_sig = [imu_sig(:, 1) blanks imu_sig(:, 2)];
                
            case 3
               
        end
        
        str_curr_temp = strsplit(errorJoints_RLDefault{j}{1}, '_');
        str_curr = [str_curr_temp{2}];
%         rmse_rldefault_tmp = rmseFct2(mocap_sig, imu_sig);
        
%         for k = 1:length(rmse_rldefault_tmp)
%             if rmse_rldefault_tmp(k) < 1e-3
%                 rmse_rldefault_tmp(k) = 0;
%             end
%         end
%         
%         indsPlot = find(rmse_rldefault_tmp > 0);
%         meanVal = sum(rmse_rldefault_tmp)/length(indsPlot);
%         
%         rmse_rldefault(j, :) = [rmse_rldefault_tmp meanVal];
%         dat_rldefault(j, :) = {str_curr ...
%             rmse_rldefault(j, 1) ...
%             rmse_rldefault(j, 2) ...
%             rmse_rldefault(j, 3) ...
%             rmse_rldefault(j, 4)};
        
%         if j > 4
%             subplotind = j+1;
%         else
            subplotind = j;
%         end
        
        ax_frame(j) = subplot(2, 4, subplotind);
        hold on;
        plot(mocap_time, mocap_sig(:,1), 'r-', 'DisplayName', 'mocap x');
        plot(mocap_time, mocap_sig(:,2), 'g-', 'DisplayName', 'mocap y');
        plot(mocap_time, mocap_sig(:,3), 'b-', 'DisplayName', 'mocap z');
%         plot(imu_time, imu_sig(:,1), 'r--', 'DisplayName', 'imu x');
%         plot(imu_time, imu_sig(:,2), 'g--', 'DisplayName', 'imu y');
%         plot(imu_time, imu_sig(:,3), 'b--', 'DisplayName', 'imu z');
        xlabel('[s]');
        
        switch subplotind
            case 1
                ylabel('[deg]');
            otherwise
                ylabel('[mm?]');
        end
        
        title([str_curr]);
    end
    legend('show');
    linkaxes(ax_frame, 'x');
end