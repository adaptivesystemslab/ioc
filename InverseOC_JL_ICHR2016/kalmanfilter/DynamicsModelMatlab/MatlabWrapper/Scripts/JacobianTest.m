%Just some stuff to test Jacobians

%Add main folder
addpath('..\');
mdlFilepath = '..\Models\avatar_v04_p_bothLegs_zyx_hip_no_ank.xml';
%Create Model and visualize
mdl = rlCModel(mdlFilepath);
mdl.forwardPosition
sens = SensorCore('sens');
sens.addDecorator('accelerometer');
mdl.addSensor(sens,'post:b_left_foot',eye(4));
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

%mdl.base = 'post:b_right_foot';
mdl.position = rand(numel(mdl.velocity));
%mdl.velocity = rand(numel(mdl.velocity))*10;
mdl.acceleration = rand(numel(mdl.velocity))*10;

mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

x = [mdl.position;mdl.velocity;mdl.acceleration];
dof = numel(mdl.joints);
%Get predicted measurement given current state
z = sens.measurement;
H_num = zeros(numel(z),numel(x));
%We modify the state a bit at a time and get new measurement
eps = 1e-8;
for i=1:numel(x)
    %Change state a little bit
    x(i) = x(i) + eps;
    %Get new measurement
    mdl.position = x(1:dof);
    mdl.velocity = x(dof+1:dof*2);
    mdl.acceleration = x(dof*2+1:end);
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    z_new = sens.measurement;
    H_num(:,i) = (z_new-z)/eps;
    %Change state back
    x(i) = x(i)-eps;
end

mdl.position = x(1:dof);
mdl.velocity = x(dof+1:dof*2);
mdl.acceleration = x(dof*2+1:end);
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.calculateSensorJacobians();
H = sens.obsJacobian;

