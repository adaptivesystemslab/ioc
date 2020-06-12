%This will reade the XML model file and load a rlCModel Object
model = rlCModel('Models\Lower_Body_Torso.xml');

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
M2 = SensorCore('ElbowSens');
M2.addDecorator('accelerometer');
M2.addDecorator('gyroscope');


%Transformation matrix defining rigid attachement of the sensor (10cm in x)
T = eye(4); T(1:2,4) = 0.1;
model.addSensor(M1,'rknee0',T);
model.addSensor(M2,'rankle0',T);

%Run forward kinematics (position only). This is needed for correct
%visualization of sensors and model before we add it to the visualizer
model.forwardPosition();
model.base = 'rankle3';

%We will now create a visualizer window and add the model to it 
vis = rlVisualizer('VIS',640,480);
%Add the model to the scene
vis.addModel(model);
%Update visualizer
vis.update();

%10 seconds of waving sampling at 100Hz
%Our timestep is 0.01
dt = 0.01;
t = [0:0.01:10]';
%Here we simply will make all revolute joints of the model follow a cosine
q = 2*cos(2*t)-2*cos(2*t(1));
dq = -4*sin(2*t);
ddq = -8*cos(2*t);

%Set torso to zero then set hip, knee, ankle to actual values, then zero
%for all else
Q = [zeros(size(q,1),3) q zeros(size(q,1),2) q q zeros(size(q,1),7)];
DQ = [zeros(size(q,1),3) dq zeros(size(q,1),2) dq dq zeros(size(q,1),7)];;
DDQ = [zeros(size(q,1),3) ddq zeros(size(q,1),2) ddq ddq zeros(size(q,1),7)];;

%Preallocate our generated mesurements [x y z] for 3 markers
mes = zeros(numel(t),numel(vertcat(model.sensors.measurement)));

%model.base = 'rankle3';
%model.forwardPosition();


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
    mes(i,:) = vertcat(model.sensors.measurement);
    
    %pause(dt);
    %vis.update(); 
end

%First Lets reset the model so we are not cheating and starting too close
%to true value. 
model.position(:) = Q(1,:);
model.velocity(:) = 0;
model.acceleration(:) = 0;
model.forwardPosition();
model.forwardVelocity();
model.forwardAcceleration();

%Create the EKF estimator to estimate joint positions, velocities, and
%accelerations based on our arm model
model.base = 'rankle3';
model.forwardPosition();

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
eta = 0.1;
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

for i=1:numel(t)
    
   %This is the measurement for current timestep, must be vertical vector
   z = mes(i,:)';
   tic;
   if i==51
      disp('bp') 
   end
   ekf.run_iteration(dt,z);

   iteration_time(i) = toc;
   %Save the estimation
   state_est(i,:) = ekf.state;
   vis.update(); 
end


