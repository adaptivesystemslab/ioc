%This will reade the XML model file and load a rlCModel Object
model = rlCModel('Models/arm_fixed.xml');

%Next we will attach 1 marker to the shoulder before the 3 revolute
%joitns, 2 markers to the elbow right before the joint, and 1 marker to the
%wrist. 

%We create Sensor Objects giving them unique names and letting them know
%that they will be measuring position and velocity. Notice we need 2
%markers at the elbow for full observability. If we only use one the EKF
%may snap to incorrect solution
M1 = SensorCore('ShoulderSens');
M1.addDecorator('accelerometer');
M1.addDecorator('gyroscope');
M1.addDecorator('yaw');
M2 = SensorCore('ElbowSens');
M2.addDecorator('accelerometer');
M2.addDecorator('gyroscope');


%Transformation matrix defining rigid attachement of the sensor (10cm in x)
T = eye(4); T(1:2,4) = 0.1;
model.addSensor(M1,'post:rhumerus',T);
model.addSensor(M2,'post:rradius',T);

%Run forward kinematics (position only). This is needed for correct
%visualization of sensors and model before we add it to the visualizer
model.forwardPosition();

%We will now create a visualizer window and add the model to it 
%vis = rlVisualizer('VIS',640,480);
%Add the model to the scene
%vis.addModel(model);
%Update visualizer
%vis.update();

%10 seconds of waving sampling at 100Hz
%Our timestep is 0.01
dt = 0.01;
t = [0:0.01:10]';
%Here we simply will make all revolute joints of the model follow a cosine
q = 0.5*cos(2*t);
dq = -sin(2*t);
ddq = -2*cos(2*t);
%Recall our model has 3 prismatic joints to position the arm in 3d space,
%then 3 revolute joints for the shoulder and a single joint at the elbow.
%Lets position the arm at [0.5 0.5 0.5] and assume that doesnt change. So
%the combined Q is as follows
Q = [q q q q];
DQ = [dq dq dq dq];
DDQ = [ddq ddq ddq ddq];

%Preallocate our generated mesurements [x y z] for 3 markers
mes = SensorMeasurement(numel(t),numel(model.sensors));

%Lets use this Q to wave the arm around, visualize, and generate
%measurements.

for i=1:numel(t)
   
    %Set the joint angles
    model.position(:) = Q(i,:)';
    %Set joint velocities
    model.velocity(:) = DQ(i,:)';
    model.acceleration(:) = DDQ(i,:)';
    
    %Run forward position and velocity, we do not run forwardAcceleration since our markers only measure position and velocity. 
    model.forwardPosition();
	model.forwardVelocity();
    model.forwardAcceleration();
    
    %Grab the generated sensor measurements.
    mes(i,:) = SensorMeasurement(model.sensors);
    
    %pause(dt);
    %vis.update(); 
end

%First Lets reset the model
model.position = Q(1,:)';
model.velocity(:) = DQ(1,:);
model.acceleration(:) = DDQ(1,:);
model.forwardPosition();
model.forwardVelocity();
model.forwardAcceleration();

%Create the EKF estimator to estimate joint positions, velocities, and
%accelerations based on our arm model
ekf = EKF_Q_DQ_DDQ(model);

%Preallocate mem for estimated [joint angles, velocities, acclerations]
state_est = zeros(numel(t),numel(ekf.state));

%Preallocate mem for timing
iteration_time = zeros(numel(t),1);

%EKF noise parameters that require tuning for best performance

%Measusrement noise, set to very low says our markers are very accurate
%Originally the noise is set to identity so here just decrease it
ekf.observation_noise = ekf.observation_noise*0.01;

%Process noise of EKF left at identity
%This eta is result of the optimization
dim = numel(model.position);
eta = 1;
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
ekf.process_noise = P;

matches = [1:numel(model.sensors)]';
matches = [matches matches];

for i=1:numel(t)
    
   %This is the measurement for current timestep, must be vertical vector
   z = mes(i,:)';
   tic;

   ekf.run_iteration(dt,z,matches);

   iteration_time(i) = toc;
   %Save the estimation
   state_est(i,:) = ekf.state;
   %vis.update(); 
end

disp('Done');
plot(state_est(:,1:2));
hold on;
plot(q,'--');
