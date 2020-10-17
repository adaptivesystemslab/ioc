function cost = half_joint_task_objective(mdl, q, dq, ddq, index1, index2)
% Purpose: compute the cost for the half_joint_half_task_objective
% Input: 2 time indices
% Prereq: run process_data.m so that the following are loaded
% - mdl, the RL model
% - POSITION - joint positions of the trajectory

% Apply the joint angles for 1st time index
[pos1, q1] = getXandQ(mdl, q, dq, ddq, index1);
[pos2, q2] = getXandQ(mdl, q, dq, ddq, index2);

% compute the displacement in task space
displace = (pos1 - pos2);

% compute the distance in joint space
distance = sqrt((q1(1)-q2(1))^2 + (q1(2)-q2(2))^2 + (q1(3)-q2(3))^2);

% Compute the cost for the half_joint_half_task_objective
cost = distance * 0.5 + norm(displace) * 0.5;
