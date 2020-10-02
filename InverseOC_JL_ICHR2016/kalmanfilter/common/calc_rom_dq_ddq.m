function [jointRom, jointMaxDq, jointMaxDdq] = calc_rom_dq_ddq(featureSet)
%     for ind_joints = 1:size(featureSet.q, 2)
%         jointRom(ind_joints) = range(featureSet.q(:, ind_joints));
%         jointMaxDq(ind_joints) = max(featureSet.dq(:, ind_joints));
%         jointMaxDdq(ind_joints) = max(featureSet.ddq(:, ind_joints));
%     end

    jointRom = range(featureSet.q);
    jointMaxDq = max(abs(featureSet.dq));
    jointMaxDdq = max(abs(featureSet.ddq));
    
    jointRom = rad2deg(jointRom);
    jointMaxDq = rad2deg(jointMaxDq);
    jointMaxDdq = rad2deg(jointMaxDdq);
end