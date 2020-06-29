% Constrained EKF test

%Add main folder
addpath('..\');
%Create Model and visualize
mdl1 = rlCModel('simple_3dof_model.xml');

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

%% Generate Trajectory in ankle knee and hip

dt = 0.01;
t = 0:dt:20;

q1 = sin(2*t);
dq1 = 2*cos(2*t);
ddq1 = -2*2*sin(2*t);

q2 = -q1;
dq2 = -dq1;
ddq2 = -ddq1;

%Run the trajectory and get measurements

mes1 = SensorMeasurement(numel(t),numel(mdl1.sensors));
mes2 = SensorMeasurement(numel(t),numel(mdl1.sensors));


for i=1:numel(t)
    mdl1.position(:) = q1(i);
    mdl1.velocity(:) = dq1(i);
    mdl1.acceleration(:) = ddq1(i);
    mdl1.forwardPosition();
    mdl1.forwardVelocity();
    mdl1.forwardAcceleration();
    
    mdl2.position(:) = q2(i);
    mdl2.velocity(:) = dq2(i);
    mdl2.acceleration(:) = ddq2(i);
    mdl2.forwardPosition();
    mdl2.forwardVelocity();
    mdl2.forwardAcceleration();
    
    mes1(i,:) = SensorMeasurement(mdl1.sensors);
    mes2(i,:) = SensorMeasurement(mdl2.sensors);
    
    vis.update();
end


%% Try estimating using Jacobian Constrained EKF

%Reset Model 
mdl1.position(:) = q1(1);
mdl1.velocity(:) = dq1(1);
mdl1.acceleration(:) = ddq1(1);
mdl1.forwardPosition;mdl1.forwardVelocity;mdl1.forwardAcceleration;
mdl2.position(:) = q2(1);
mdl2.velocity(:) = dq2(1);
mdl2.acceleration(:) = ddq2(1);
mdl2.forwardPosition;mdl2.forwardVelocity;mdl2.forwardAcceleration;
vis.update();

ekf = DC_GP_EKF_Q_DQ_DDQ(mdl1,mdl2);

%Set up obs noise
ekf.observation_noise = diag(repmat([0.01 0.01 0.01 1 1 1],1,6));

G = [];
dim = numel(mdl1.joints);
for i=ekf.sizeX/dim/2:-1:1
    G = [G; ones(dim,1)*dt^i / factorial(i)];
end
%Set up proc noise
P_tmp = G*G'*10;
P = zeros(size(P_tmp));
for i=1:ekf.sizeX/dim/2
    for j=i:ekf.sizeX/dim/2
        P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
            diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
    end
end

%What if there is no noise in position 
%P(1:3,:) = 0;
P = P+P' - diag(diag(P));

ekf.process_noise(1:end/2,1:end/2) = P;
ekf.process_noise(end/2+1:end,end/2+1:end) = P;

ekf.covariance(:,:) = 0;
%ekf.covariance(10:12,:) = 0;

matches = [1:6]';
matches = [matches matches];
%% Do the estimation 

f_ee1 = mdl1.getFrameByName('frame_ee');
f_ee2 = mdl2.getFrameByName('frame_ee');

vel1 = zeros(numel(t),6);
vel2 = zeros(numel(t),6);
accel1 = zeros(numel(t),3);
accel2 = zeros(numel(t),3);
pos1 = zeros(numel(t),3);
rpy1 = zeros(numel(t),3);
pos2 = zeros(numel(t),3);
rpy2 = zeros(numel(t),3);
ekf_state = zeros(numel(t),ekf.sizeX);
ekf_cov_trace = zeros(numel(t),1);
g = mdl1.g;
for i=1:numel(t)
    
    mes2_arr = mes2(i,:).getMesArray;
    %if(i > 100)
        mes2_arr([3 9]) = mes2_arr([3 9]) + 1;
    %end
    ekf.run_iteration(dt,[mes1(i,:).getMesArray mes2_arr]')
    ekf_state(i,:) = ekf.state;
    vel1(i,:) = [f_ee1.t(1:3,1:3)*f_ee1.v(1:3); f_ee1.t(1:3,1:3)*f_ee1.v(4:6)]; 
    vel2(i,:) = [f_ee2.t(1:3,1:3)*f_ee2.v(1:3); f_ee2.t(1:3,1:3)*f_ee2.v(4:6)];
    accel1(i,:) = f_ee1.t(1:3,1:3)*(f_ee1.a(4:6) + cross(f_ee1.v(1:3),f_ee1.v(4:6)))-g;
    accel2(i,:) = f_ee2.t(1:3,1:3)*(f_ee2.a(4:6) + cross(f_ee2.v(1:3),f_ee2.v(4:6)))-g;
    pos1(i,:) = f_ee1.t(1:3,4);
    pos2(i,:) = f_ee2.t(1:3,4);
    ekf_cov_trace(i) = trace(blkdiag(ekf.covariance(1:3,1:3),ekf.covariance(10:12,10:12)));

    T1 = f_ee1.t;
    t1 = T1(1:3,4);
    [yaw, pitch, roll] = dcm2angle(T1(1:3,1:3)');
    r1 = [yaw pitch roll]';
    rpy1(i,:) = r1;
    T2 = f_ee2.t;
    t2 = T2(1:3,4);
    [yaw, pitch, roll] = dcm2angle(T2(1:3,1:3)');
    r2 = [yaw pitch roll]';
    rpy2(i,:) = r2;
    
    
    vis.update();
end



