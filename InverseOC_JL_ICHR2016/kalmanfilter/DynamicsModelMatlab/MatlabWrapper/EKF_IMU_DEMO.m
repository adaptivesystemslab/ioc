%DEMO OF HOW TO USE MAGIC EKF WITH IMUs 

%% Here we load up the model, attach 3 imus to it and wave it around
vis = rlVisualizer('vis',640,480);

mdl = rlCModel('Models\simple_3dof_model.xml');
mdl.forwardPosition();

T = eye(4,4);
T(1:3,4) = 0.1;

%Here we add 3 IMUs
imu1 = SensorCore('imu1');
imu1.addDecorator('gyroscope');
imu1.addDecorator('accelerometer');
mdl.addSensor(imu1,'frame1',T);

imu2 = SensorCore('imu2');
imu2.addDecorator('gyroscope');
imu2.addDecorator('accelerometer');
mdl.addSensor(imu2,'frame3',T);

imu3 = SensorCore('imu3');
imu3.addDecorator('gyroscope');
imu3.addDecorator('accelerometer');
mdl.addSensor(imu3,'frame_ee',T);
mdl.forwardPosition();

vis.addModel(mdl);
vis.update();

dt = 0.01;
t = 0:dt:10;
multiplier = 2;
q = cos(multiplier*t);
dq = -multiplier*sin(multiplier*t);
ddq = -multiplier^2*cos(multiplier*t);

%Preallocate measurement object array 
mes = SensorMeasurement(numel(t),numel(mdl.sensors));

for i=1:numel(t)
    %Set the position, velocity, and acceleration for the model
   mdl.position(:) = q(i);
   mdl.velocity(:) = dq(i);
   mdl.acceleration(:) = ddq(i);
   
   %Forward Kinematics
   mdl.forwardPosition();
   mdl.forwardVelocity();
   mdl.forwardAcceleration();
   
   %Store Sensor Measurements
    mes(i,:) = SensorMeasurement(mdl.sensors);
    vis.update();
end

%% Now lets estimate using the generated measurements

%Rest model
mdl.position(:) = q(1);
mdl.velocity(:) = dq(1);
mdl.acceleration(:) = ddq(1);

%Create EKF Object 
ekf = EKF_Q_DQ_DDQ(mdl);

%Set EKF Observation Noise properties
ACCEL_NOISE = [10 10 10];
GYRO_NOISE = [0.01 0.01 0.01];
ekf.observation_noise = diag(repmat([GYRO_NOISE ACCEL_NOISE],1,numel(mdl.sensors)));

%Set EKF Process Noise
eta = 10;
dim = numel(mdl.joints);
G = [ones(dim,1)*dt^2/2*eta; ones(dim,1)*dt*eta; ones(dim,1)*eta];
P_tmp = G*G';
P = zeros(size(P_tmp));
for i=1:3
    for j=i:3
        P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
            diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
    end
end
P = P+P' - diag(diag(P));
ekf_eul.process_noise = P;

%Tell EKF that measurements match one to one 
matches = [1 1; 2 2; 3 3];

%pre-allocate state estimate results
ekf_state_est = zeros(numel(t),ekf.sizeX);
mes_est = SensorMeasurement(numel(t),numel(mdl.sensors));
%Run EKF 
for i=1:numel(t)
    z = mes(i,:);
    ekf.run_iteration(dt,z,matches);
    ekf_state_est(i,:) = ekf.state;
    mes_est(i,:) = ekf.makeMeasure(ekf.state);
    vis.update();
end

%Convert objects to array
mes_array = mes.getMesArray();
mes_est_array = mes_est.getMesArray();