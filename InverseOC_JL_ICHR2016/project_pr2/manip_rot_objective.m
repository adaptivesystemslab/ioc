function cost = manip_rot_objective(mdl, q, dq, ddq, index1)
% Purpose: compute the cost for the manip_rot_objective
% Input: a time index
% Prereq: run process_data.m so that the following are loaded
% - mdl, the RL model
% - POSITION - joint positions of the trajectory
% - VELOCITY - joint velocities of the trajectory
% - ACCELERATION - joint accelerations of the trajectory

% % Apply the joint positions, velocties, and accelerations for the time index
mdl = inputData(mdl, q(:, index1), dq(:, index1), ddq(:, index1));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

mdl.calculateJacobian();

% Grab the lower 3 rows of the Jacobian (rotational component)
J2 = mdl.J(4:6,:); 

% Compute the cost (mdl.J is the Jacobian)
cost = 1 / sqrt(det(J2 * J2'));