%Test the multi model EKF 
%Add main folder
addpath('..\');

%Save basic video
save_video = true;

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

mdl2 = rlCModel('simple_4dof_model.xml');

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
mdl2.addSensor(sens,'frame4',eye(4));
sens = SensorCore('sens3');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');

%Align last sensors
mdl2.addSensor(sens,'frame6',eye(4));

%Place second model so last sensors align
mdl2.transforms(1).t(1:3,1:3) = rotz(180)*mdl2.transforms(1).t(1:3,1:3);
mdl2.transforms(1).t(1:3,4) = [0 -0.4 0]';
mdl2.forwardPosition();


vis = rlVisualizer('vis',1280,960);
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
    
    mdl2.position(2:end) = q2(i);
    mdl2.velocity(2:end) = dq2(i);
    mdl2.acceleration(2:end) = ddq2(i);
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
mdl2.position(2:end) = q2(1);
mdl2.velocity(2:end) = dq2(1);
mdl2.acceleration(2:end) = ddq2(1);
mdl2.position(1)=0;mdl2.velocity(1)=0;mdl2.acceleration(1)=0;
mdl2.forwardPosition;mdl2.forwardVelocity;mdl2.forwardAcceleration;
vis.update();

ekf = MC_EKF_Q_DQ_DDQ();
%Add the models
ekf.addModel(mdl1);
ekf.addModel(mdl2);
%Set up obs noise
ekf.observation_noise = diag(repmat([0.1 0.1 0.1 0.1 0.1 0.1],1,6));
%ekf.covariance(:,:) = 0;
%Add constraints that the end effectors are in the same position 
%C1 = EKF_Constraint_BF('C1',mdl1,'frame_ee',eye(4),mdl2,'frame_ee',eye(4));
%C1.type(:) = false;C1.type([1 2 3 4 5 6]) = true;
%C2 = EKF_Constraint_T0EE('C2',mdl1,'frame_ee',eye(4),mdl2,'frame_ee',eye(4));
%C2.type(:) = false;C2.type([1]) = true;
mdl1.forwardPosition();
mdl2.forwardPosition();
%ekf.addConstraint(C1);
%ekf.addConstraint(C2);

%Because constraint attaches more sensors we restart visualizer
vis = rlVisualizer('vis',1280,960);
vis.addModel(mdl1);
vis.addModel(mdl2);
vis.update();
vis.setBackgroundColour([1 1 1])
vis.setViewportPose([-0.4828   -1.0186         0   53.0000   15.6000    3.2859]');
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

mes_arr = [mes1.getMesArray mes2.getMesArray];
mes_noise = randn(size(mes_arr)).*repmat([0.1 0.1 0.1 1 1 1],size(mes_arr,1),6);
mes_arr_noisy = mes_arr + mes_noise;
mes_arr_noisy_drift = mes_arr_noisy;
%mes_arr_noisy_drift(:,19:21) = mes_arr_noisy_drift(:,19:21)+repmat(randn(1,3)*0.5,size(mes_arr,1),1);
%mes_arr_noisy_drift(:,25:27) = mes_arr_noisy_drift(:,25:27)+repmat(randn(1,3)*0.5,size(mes_arr,1),1);
mes_arr_noisy_drift(:,[1:3 19:21 25:27 31:33]) = mes_arr_noisy_drift(:,[1:3 19:21 25:27 31:33])+0.1;
mes_arr_noisy_drift(:,[20]) = mes_arr_noisy_drift(:,[20])+1;

if save_video
   video_path = 'C:\asl_git\kalmanfilter\General_FKEKF\DynamicsModelMatlab\MatlabWrapper\results\BF_C.avi';
   v = VideoWriter(video_path); 
   v.open();
end

g = mdl1.g;
for i=1:numel(t)
    ekf.run_iteration(dt,mes_arr_noisy_drift(i,:)')
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
    if save_video
       [img,imgbool] = vis.getScreenshot();
       if imgbool
           writeVideo(v,img);
       end
    end
end

if save_video
   v.close(); 
end

%% Test the axis alignment constraint 

clear vis;
clear mdl;

mdl1 = rlCModel('..\Models\spherical.xml');
mdl1.forwardPosition();

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

mdl2 = rlCModel('..\Models\spherical2.xml');

%add sensor to end effector
sens = SensorCore('sens1');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');
%sens.addDecorator('velocity');
T = eye(4);
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

q2 = q1;
dq2 = dq1;
ddq2 = ddq1;

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

ekf = MC_EKF_Q_DQ_DDQ();
%Add the models
ekf.addModel(mdl1);
ekf.addModel(mdl2);
%Set up obs noise
ekf.observation_noise = diag(repmat([0.01 0.01 0.01 1 1 1],1,6));
%ekf.covariance(:,:) = 0;
%Add constraints that the end effectors are in the same position 
C = EKF_Constraint_T0EE('C1',mdl1,'frame_ee',eye(4),mdl2,'frame_ee',eye(4));
C.type = logical([0 0 1 0 0 0]');
mdl1.forwardPosition();
mdl2.forwardPosition();
ekf.addConstraint(C);

%Because constraint attaches more sensors we restart visualizer
vis = rlVisualizer('vis',1280,960);
vis.addModel(mdl1);
vis.addModel(mdl2);
vis.setBackgroundColour([1 1 1]');
vis.setViewportPose([-0.1867   -0.2009         0   29.0000   29.0000    3.3339]');
vis.update();

%% Do the estimation 
mdl1.position(:) = q1(1);
mdl1.velocity(:) = dq1(1);
mdl1.acceleration(:) = ddq1(1);
mdl1.forwardPosition;mdl1.forwardVelocity;mdl1.forwardAcceleration;
mdl2.position(:) = q2(1);
mdl2.velocity(:) = dq2(1);
mdl2.acceleration(:) = ddq2(1);
mdl2.forwardPosition;mdl2.forwardVelocity;mdl2.forwardAcceleration;
vis.update();

if save_video
   video_path = 'C:\asl_git\kalmanfilter\General_FKEKF\DynamicsModelMatlab\MatlabWrapper\results\EE_C.avi';
   v = VideoWriter(video_path); 
   v.open(); 
end


f_ee1 = mdl1.getFrameByName('frame_ee');
f_ee2 = mdl2.getFrameByName('frame_ee');

vel1 = zeros(numel(t),6);
vel2 = zeros(numel(t),6);
accel1 = zeros(numel(t),3);
accel2 = zeros(numel(t),3);
pos1 = zeros(numel(t),3);
rpy1 = zeros(numel(t),3);
zax1 = zeros(numel(t),3);
pos2 = zeros(numel(t),3);
rpy2 = zeros(numel(t),3);
zax2 = zeros(numel(t),3);
t12 = zeros(numel(t),3);
ekf_state = zeros(numel(t),ekf.sizeX);
ekf_cov_trace = zeros(numel(t),1);
g = mdl1.g;
for i=1:numel(t)
    
    mes2_arr = mes2(i,:).getMesArray;
    %if(i > 100)
        mes2_arr([1 3 9]) = mes2_arr([1 3 9]) + 1;
        mes2_arr([4 10]) = mes2_arr([4 10]) + 0.3;
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
    zax1(i,:) = T1(1:3,3);
    t1 = T1(1:3,4);
    [yaw, pitch, roll] = dcm2angle(T1(1:3,1:3)');
    r1 = [yaw pitch roll]';
    rpy1(i,:) = r1;
    T2 = f_ee2.t;
    zax2(i,:) = T2(1:3,3);
    t2 = T2(1:3,4);
    [yaw, pitch, roll] = dcm2angle(T2(1:3,1:3)');
    r2 = [yaw pitch roll]';
    rpy2(i,:) = r2;
    
    T12 = SE3.fastInverse(T1)*T2;
    t12(i,:) = T12(1:3,4);
    vis.update();
    if save_video
       [img,imgbool] = vis.getScreenshot();
       if imgbool
           writeVideo(v,img);
       end
    end
end

if save_video
   v.close(); 
end

