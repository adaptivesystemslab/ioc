function [q, dq, ddq, dddq] = calc_joint_trajectory(spline_fit, x_traj)
% using splines, calculate the components needed for the 
% x_traj_to_use = [x_traj-h; x_traj; x_traj+h];
% x_traj_to_use = x_traj_to_use(:, x_knot_ind); % remove the edge pieces to prevent extrapolation effects
x_traj_to_use = x_traj;

q = [];
dq = [];
ddq = [];
dddq = [];
for ii = 1:size(spline_fit, 1) %number of joints
    q_slice =    fnval(spline_fit(ii).sp_qd0,x_traj_to_use);
    dq_slice =   fnval(spline_fit(ii).sp_qd1,x_traj_to_use);
    ddq_slice =  fnval(spline_fit(ii).sp_qd2,x_traj_to_use);
    dddq_slice = fnval(spline_fit(ii).sp_qd3,x_traj_to_use);
    
    q(:, :, ii)    = q_slice;
    dq(:, :, ii)   = dq_slice;
    ddq(:, :, ii)  = ddq_slice;
    dddq(:, :, ii) = dddq_slice;
%     
%     q_m(ii, :, :)    = reshape(fnval(spline_fit(ii).sp_qd0,x_traj_to_use), [size(x_traj_to_use, 2) size(x_traj_to_use, 1)]);
%     dq_m(ii, :, :)   = fnval(spline_fit(ii).sp_qd1,x_traj_to_use);
%     ddq_m(ii, :, :)  = fnval(spline_fit(ii).sp_qd2,x_traj_to_use);
%     dddq_m(ii, :, :) = fnval(spline_fit(ii).sp_qd3,x_traj_to_use);
end
end