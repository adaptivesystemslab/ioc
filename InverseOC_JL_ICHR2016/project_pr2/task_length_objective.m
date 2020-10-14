function cost = task_length_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the task_length_objective
% Input: 2 time indices
% Prereq: run process_data.m (mdl and POSITION need to be loaded)

% Apply the joint angles for 1st time index
mdl = inputData(mdl, q(:, index1), dq(:, index1), ddq(:, index1));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

% Store the homogenous transform of the end effector (frame 13, gripper_tool_frame)
end_eff = mdl.getFrameByName('frame13'); 

% Get the position of the end effector
pos1 = end_eff.t(1:3, 4);

% Repeat the same process with 2nd time index
mdl = inputData(mdl, q(:, index2), dq(:, index2), ddq(:, index2));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

end_eff = mdl.getFrameByName('frame13');
pos2 = end_eff.t(1:3, 4);

% Compute the cost for the task_length_objective
displace = (pos1 - pos2);
cost = norm(displace);