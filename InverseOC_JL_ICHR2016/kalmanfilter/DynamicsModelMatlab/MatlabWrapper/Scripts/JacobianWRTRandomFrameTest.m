%Test Jacobians WRT random frame

addpath('..\');
%Create Model and visualize
mdl = rlCModel('..\Models\Lower_Body_Rev.xml');
mdl.forwardPosition();

%Add Sensor to the knee
sens = SensorCore('knee_imu');
sens.addDecorator('yaw');

knee_frame = mdl.getFrameByName('rknee0');
asis_frame = mdl.getFrameByName('mid_asis');
T = eye(4);
T(1:3,1:3) = knee_frame.t(1:3,1:3)';
T(1:3,4) = T(1:3,1:3)*[0.1 0 0.1]';
mdl.addSensor(sens,'rknee0',T);

mdl.forwardPosition();

%Set the sensor's base to be the hip
sens.base = asis_frame;

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();


