function cost = task_jerk_objective(mdl, q, dq, ddq)
% Purpose: compute the cost for the task_jerk
% Input: 1 time index

% Apply the joint angles for the time index
for i = 1:size(q, 2)
    [pos, ~, ~] = getXandQ(mdl, q, dq, ddq, i);
    x(:,i) = pos;
end

dt = 1/100;
x =  filter_dualpassBW(x', 0.05, 0, 4)';

dx = calcDeriv(x, dt);
ddx = calcDeriv(dx, dt);
dddx = calcDeriv(ddx, dt);

cost = dddx;