function cost = orientation_length_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the orientation_length_objective
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

% Get the rotation of the end effector
rot1 = end_eff.t(1:3,1:3);

% Repeat the same process with 2nd time index
mdl = inputData(mdl, q(:, index2), dq(:, index2), ddq(:, index2));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

end_eff = mdl.getFrameByName('frame13');
rot2 = end_eff.t(1:3,1:3);

% Compute the cost for the orientation_length_objective
arg = (trace(rot1 *  rot2') - 1) / 2;

% force a round if it's above 1 or below -1 to prevent complex values
% if arg > 1
%     arg = 1;
% end
% 
% if arg < -1
%     arg = -1;
% end

cost = acos(arg);
