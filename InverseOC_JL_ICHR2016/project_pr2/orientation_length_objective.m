function cost = orientation_length_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the orientation_length_objective
% Input: 2 time indices
% Prereq: run process_data.m so that the following are loaded
% - mdl, the RL model
% - POSITION - joint positions of the trajectory

% Apply the joint angles for 1st time index
[~, ~, rot1] = getXandQ(mdl, q, dq, ddq, index1);
[~, ~, rot2] = getXandQ(mdl, q, dq, ddq, index2);

% Compute the cost for the orientation_length_objective
arg = (trace(rot1 *  rot2') - 1) / 2;

cost = acos(arg);
