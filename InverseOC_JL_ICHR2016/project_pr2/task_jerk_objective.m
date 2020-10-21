function cost = task_jerk_objective(mdl, q, dq, ddq, index)
% Purpose: compute the cost for the task_jerk
% Input: 1 time index
% Prereq: run process_data.m (mdl and POSITION need to be loaded)

% Apply the joint angles for the time index
[pos, ~, ~] = getXandQ(mdl, q, dq, ddq, index);
x = pos;
dt = 1/100;

dx = calcDeriv(x, dt);
ddx = calcDeriv(dx, dt);
dddx = calcDeriv(ddx, dt);

cost = dddx;