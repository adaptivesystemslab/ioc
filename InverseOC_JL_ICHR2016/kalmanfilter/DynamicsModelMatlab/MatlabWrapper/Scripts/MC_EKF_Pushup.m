addpath('..\');
%Create Model and visualize
mdl1 = rlCModel('..\Models\planar2.xml');

%add sensor to end effector
sens = SensorCore('sens1');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');
mdl1.addSensor(sens,'frame_ee',eye(4));
sens = SensorCore('sens2');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');
mdl1.addSensor(sens,'frame2',eye(4));
sens = SensorCore('sens3');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');
mdl1.addSensor(sens,'frame4',eye(4));
mdl1.forwardPosition();

mdl2 = rlCModel('simple_3dof_model_2.xml');

%add sensor to end effector
sens = SensorCore('sens1');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');
T = eye(4);
T(1:3,1:3) = roty(180);
mdl2.addSensor(sens,'frame_ee',T);
sens = SensorCore('sens2');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');
mdl2.addSensor(sens,'frame2',eye(4));
sens = SensorCore('sens3');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');

%Align last sensors
mdl2.addSensor(sens,'frame4',eye(4));

%Place second model so last sensors align
mdl2.transforms(1).t(1:3,1:3) = rotz(180)*mdl2.transforms(1).t(1:3,1:3);
mdl2.transforms(1).t(1:3,4) = [0 -0.4 0]';
mdl2.forwardPosition();


vis = rlVisualizer('vis',640,480);
vis.addModel(mdl1);
vis.addModel(mdl2);
vis.update();
