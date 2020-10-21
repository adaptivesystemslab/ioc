function cost = joint_jerk_objective(mdl, q, dq, ddq, index)
% Purpose: compute the cost for the joint_jerk
% Input: 1 time index
% Prereq: run process_data.m (mdl and POSITION need to be loaded)

[~, qi] = getXandQ(mdl, q, dq, ddq, index);
dt = 1/100;

dqi = calcDeriv(qi, dt);
ddqi = calcDeriv(dqi, dt);
dddqi = calcDeriv(ddqi, dt);

cost = dddqi;