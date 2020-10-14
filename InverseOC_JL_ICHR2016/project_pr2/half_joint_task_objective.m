function cost = half_joint_task_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the half_joint_half_task_objective
% Input: 2 time indices
% Prereq: run process_data.m so that the following are loaded
% - mdl, the RL model
% - POSITION - joint positions of the trajectory

% Apply the joint angles for 1st time index
mdl = inputData(mdl, q(:, index1), dq(:, index1), ddq(:, index1));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

% Store the homogenous transform of the end effector (frame 13, gripper_tool_frame)
end_eff = mdl.getFrameByName('frame13'); 
q1 = q(:, index1);
pos1 = end_eff.t(1:3, 4);

% Repeat the same process with 2nd time index
mdl = inputData(mdl, q(:, index2), dq(:, index2), ddq(:, index2));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

end_eff = mdl.getFrameByName('frame13');
q2 = q(:, index2);
pos2 = end_eff.t(1:3, 4);

% compute the displacement in task space
displace = (pos1 - pos2);

% compute hte distance in joint space
distance = sqrt((q1(1)-q2(1))^2 + (q1(2)-q2(2))^2 + (q1(3)-q2(3))^2);

% Compute the cost for the half_joint_half_task_objective
cost = distance * 0.5 + norm(displace) * 0.5;
