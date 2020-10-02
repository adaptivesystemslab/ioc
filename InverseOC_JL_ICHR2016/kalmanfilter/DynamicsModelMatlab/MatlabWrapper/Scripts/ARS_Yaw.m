% This is a simple script to figure out how to attach yaw sensors to ARS
% model

addpath('../');

mdl = rlCModel('C:\Cardon\ars\core\assets\avatar_v04_p_bothLegs_rev01.xml');
mdl.forwardPosition();

sens = SensorCore('yaw');
sens.addDecorator('yaw');
T = eye(4);

Rot_f_to_w = mdl.getFrameByName('post:b_spine_3').t;
Rot_f_to_w(1:3,4) = 0;
Rot_w_to_f = SE3.fastInverse(Rot_f_to_w);
Rot_b_to_w = mdl.getFrameByName('post:b_pelvis:rotateZ').t;
Rot_b_to_w(1:3,4) = 0;
T_rx = eye(4);
T_rx(1:3,1:3) = rotx(-45,'deg');
R_f_to_sens = Rot_w_to_f*Rot_b_to_w*T_rx;
R_f_to_sens(1:3,4) = 0.2;

mdl.addSensor(sens,'post:b_spine_3',R_f_to_sens);
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();
