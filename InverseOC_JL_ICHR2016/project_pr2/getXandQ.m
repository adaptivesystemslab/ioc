function [pos1, ang1, rot1] = getXandQ(mdl, q, dq, ddq, index1)
mdl = inputData(mdl, q(:, index1), dq(:, index1), ddq(:, index1));
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();

% Store the homogenous transform of the end effector (frame 13, gripper_tool_frame)
end_eff = mdl.getFrameByName('frame13'); 
ang1 = q(:, index1);
pos1 = end_eff.t(1:3, 4);
rot1 = end_eff.t(1:3, 1:3);