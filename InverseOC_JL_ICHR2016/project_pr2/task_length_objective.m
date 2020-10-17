function cost = task_length_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the task_length_objective
% Input: 2 time indices
% Prereq: run process_data.m (mdl and POSITION need to be loaded)

% Apply the joint angles for 1st time index
[pos1, ~, ~] = getXandQ(mdl, q, dq, ddq, index1);
[pos2, ~, ~] = getXandQ(mdl, q, dq, ddq, index2);

% Compute the cost for the task_length_objective
displace = (pos1 - pos2);
cost = norm(displace.^2);