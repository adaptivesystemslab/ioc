%Test parllelized Jacobian computation 

mdl = rlCModel('Models/arm_fixed.xml');
sens = SensorCore('IMU');
sens.addDecorator('accelerometer');
mdl.forwardPosition();
mdl.addSensor(sens,'post:rradius',eye(4));
mdl.forwardPosition();
mdl.g = [0 0 0]';

time_single = 0;
for i=1:1000
    tic;
    mdl.calculateSensorJacobians();
    time_single = time_single+toc;
end
disp('One Sensor Time: ');
disp(time_single/1000);


mdl = rlCModel('Models/arm_fixed.xml');
mdl.forwardPosition();
for i=1:10
    sens = SensorCore(['IMU_' num2str(i)]);
    sens.addDecorator('accelerometer');
    mdl.addSensor(sens,'post:rradius',eye(4));
end
mdl.forwardPosition();
mdl.g = [0 0 0]';

time_multi = 0;
for i=1:1000
    tic;
    mdl.calculateSensorJacobians();
    time_multi = time_multi+toc;
end
disp('Multi Sensor Time: ');
disp(time_multi/1000);



