function h = metrics_mocap_IMU(outputStruct_mocap, outputStruct_imu)
%     startJointStr = 'frame_6dof_revz0_to_frame_6dof_revz1';
%     endJointStr = 'joint_lankle_2';
%     allJointStr = outputStruct_mocap.jointNames;
%     startJointInd = find(ismember(allJointStr, startJointStr), 1);
%     endJointInd = find(ismember(allJointStr, endJointStr), 1);
%     jointInds = startJointInd:endJointInd;
% 
%     jointsToUse = {};
%     for i = 1:length(jointInds)
%         jointsToUse{i} = allJointStr{jointInds(i)};
%     end

    jointsToUse{1} = 'joint_rhip_0';
    jointsToUse{2} = 'joint_rhip_1';
    jointsToUse{3} = 'joint_rhip_2';
    jointsToUse{4} = 'joint_rknee_0';
    jointsToUse{5} = 'joint_rknee_1';
    jointsToUse{6} = 'joint_lhip_0';
    jointsToUse{7} = 'joint_lhip_1';
    jointsToUse{8} = 'joint_lhip_2';
    jointsToUse{9} = 'joint_lknee_0';
    jointsToUse{10} = 'joint_lknee_1';
    jointsToUse{11} = 'frame_6dof_revz0_to_frame_6dof_revz1';
    jointsToUse{12} = 'frame_6dof_revy0_to_frame_6dof_revy1';
    jointsToUse{13} = 'frame_6dof_revx0_to_frame_6dof_revx1';
    
    allJointStrMocap = outputStruct_mocap.jointNames;
    allJointStrImu = outputStruct_imu.jointNames;
    
    h = figure;
    for i = 1:length(jointsToUse)
        ax(i) = subplot(3, 5, i);
        currJointStr = jointsToUse{i};
        
        currJointInd = find(ismember(allJointStrMocap, currJointStr), 1);
        plot(outputStruct_mocap.time, outputStruct_mocap.q(:, currJointInd), 'b', 'DisplayName', 'mocap'); 
        
        hold on;
        
        currJointInd = find(ismember(allJointStrImu, currJointStr), 1);
        if ~isempty(currJointInd)
            plot(outputStruct_imu.time - outputStruct_imu.time(1), outputStruct_imu.q(:, currJointInd), 'r', 'DisplayName', 'IMU');
        end
        
        ylim([-pi pi]);
        
        title(currJointStr);
    end
    
    linkaxes(ax, 'x');
end
