function cost = joint_length_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the joint_length_objective
% Input: 2 time indices
% Prereq: run process_data.m so that the following are loaded
% - mdl, the RL model
% - POSITION - joint positions of the trajectory

% Apply the joint angles for 1st time index
pos1 = q(:, index1);
pos2 = dq(:, index2);

% Compute the distance in joint space
distance = sum((pos1 - pos2).^2);

% Return the cost for the joint_length_objective
cost = distance;