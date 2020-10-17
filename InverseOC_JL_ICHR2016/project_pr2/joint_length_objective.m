function cost = joint_length_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the joint_length_objective
% Input: 2 time indices
% Prereq: run process_data.m so that the following are loaded
% - mdl, the RL model
% - POSITION - joint positions of the trajectory

% Apply the joint angles for 1st time index
[~, q1] = getXandQ(mdl, q, dq, ddq, index1);
[~, q2] = getXandQ(mdl, q, dq, ddq, index2);

% Compute the distance in joint space
distance = sum((q1 - q2).^2);

% Return the cost for the joint_length_objective
cost = distance;