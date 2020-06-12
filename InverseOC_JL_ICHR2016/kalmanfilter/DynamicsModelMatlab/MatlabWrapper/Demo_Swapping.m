%Demo of Marker Based pose estimation.
%
% First we load a single arm model that can be positioned anywhere in
% space. The arm has 3 revolute joints in the shoulder simulating a ball
% joint and a single joint at the elbow. It is connected to the world frame
% through 3 prismatic joints allowing arbitrary positioning in space. We
% can modify the model to include more joints for the hand if needed. 

%It is important to clear variables carefully or matlab will crash. Many of
%the objects are in C++ and if we don't clear Matlab's memory management
%messes up. Always clear model last. Also clear these variables before
%exiting matlab or it may crash. 
clear ekf   %ekf is matlab handle but hold the model so needs to be cleared
clear vis   %visializer must be cleared before model
clear M1 M2 M3 M4   %clear sensor handles
clear model         %clear model object last


%This will reade the XML model file and load a rlCModel Object
model = rlCModel('Models\arm.xml');

%Next we will attach 1 marker to the shoulder before the 3 revolute
%joitns, 2 markers to the elbow right before the joint, and 1 marker to the
%wrist. 

%We create Sensor Objects giving them unique names and letting them know
%that they will be measuring position and velocity. Notice we need 2
%markers at the elbow for full observability. If we only use one the EKF
%may snap to incorrect solution
M1 = SensorCore('ShoulderMarker');
M1.addDecorator('position');
M2 = SensorCore('ElbowMarker1');
M2.addDecorator('position');
M3 = SensorCore('ElbowMarker2');
M3.addDecorator('position');
M4 = SensorCore('WristMarker');
M4.addDecorator('position');

%We attach the sensors to the model, take a look at the XML model file. We
%can see that after first 3 prismatic joints we have a frame called
%"rhumerus". We will attach the shoulder marker to it with a small offset.
%Note, the offset is in the frame you are attaching the sensor to. 

%Transformation matrix defining rigid attachement of the sensor (10cm in x)
T = eye(4); T(1:2,4) = 0.1;
model.addSensor(M1,'rhumerus',T);
 
%Similarly we will add M2 and M3 to the other frames on the body 
model.addSensor(M2,'rradius',T);
%Similarly we will add M2 and M3 to the other frames on the body 
T(1:2,4) = -0.1;
model.addSensor(M3,'rradius',T);

%This is actually out wrist, the naming convention is strange. We can
%change this to make it easier. 
model.addSensor(M4,'post:rradius',T);

%Run forward kinematics (position only). This is needed for correct
%visualization of sensors and model before we add it to the visualizer
model.forwardPosition();

%We will now create a visualizer window and add the model to it 
vis = rlVisualizer('VIS',640,480);
%Add the model to the scene
vis.addModel(model);
%Update visualizer
vis.update();

disp('Added Model with markers to visualizer, press any key to continue');
pause();

%So now we have a model of the arm. Next we wave it around and generate
%marker measurements. We will later use the generated marker measurements
%to estimate the joint angles.

%10 seconds of waving sampling at 100Hz
dt = 0.01;
t = [0:dt:10]';
%Here we simply will make all revolute joints of the model follow a cosine
q = 0.5*cos(2*t);
dq = -sin(2*t);
%Recall our model has 3 prismatic joints to position the arm in 3d space,
%then 3 revolute joints for the shoulder and a single joint at the elbow.
%Lets position the arm at [0.5 0.5 0.5] and assume that doesnt change. So
%the combined Q is as follows
Q = [ones(numel(t),3)*0.5 q q q q];
DQ = [zeros(numel(t),3) dq dq dq dq];

%Preallocate our generated mesurements [x y z] for 3 markers
mes = SensorMeasurement(numel(t),numel(model.sensors));

%Lets use this Q to wave the arm around, visualize, and generate
%measurements.

for i=1:numel(t)
   
    %Set the joint angles
    model.position(:) = Q(i,:)';
    %Set joint velocities
    model.velocity(:) = DQ(i,:)';
    
    %Run forward position and velocity, we do not run forwardAcceleration since our markers only measure position and velocity. 
    model.forwardPosition();
	model.forwardVelocity();
    
    %Grab the generated marker measurements.
	%This vector will be in the order that we have added the markers and decorators 
	%| Marker 1 position|
	%| Marker 1 velocity|
	%| Marker 2 position|
	%| Marker 2 velocity|
	%|.|
	%|.|
    mes(i,:) = SensorMeasurement(model.sensors);
    
    %Make BAD MEASUREMENTS
    if i > 200 && i < 300
        mes(i,end-1).mes = [10 10 10]';
    end
    
    %MAKE SWAP
    if i > 600 && i < 700
        tmp = mes(i,end).mes;
        mes(i,end).mes = mes(i,end-1).mes;
        mes(i,end-1).mes = tmp;
    end
    
    %Visualize
    %vis.update();
    %Pause for dt
    %pause(0.01);
end

%Add some noise
mes_mat = mes.getMesArray();
mes_mat = mes_mat + 0.01*randn(size(mes_mat));
mes.setMesArray(mes_mat);

disp('Done Waving arm around press any key to conitnue');
pause();

%Finally we will go back the other way, we will use the generated
%measurements to estimate the joint position, velocity, and acceleration.
%TODO

%First Lets reset the model so we are not cheating and starting too close
%to true value. 
model.position = Q(:,1)';
model.forwardPosition();
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
ekf.observation_noise = ekf.observation_noise*0.02;

%FOR NOW MAKE NOISE THE SAME
eta = 1;

%Process noise
dim = numel(model.joints);
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

%Our timestep is 0.01
dt = 0.01;

for i=1:numel(t)
    
   %This is the measurement for current timestep, must be vertical vector
   	%| Marker 1 position|
	%| Marker 1 velocity|
	%| Marker 2 position|
	%| Marker 2 velocity|
	%|.|
	%|.|
   z = mes(i,:);
   %Here we run the itiration and time it. Once you have your model set up
   %and you create the EKF_Q_DQ_DDQ object you will have to call
   %ekf.run_iteration(dt, measurement_vector) as the marker data comes in
   %from motion capture.
   tic;
   
   ekf.run_iteration(dt,z);
   iteration_time(i) = toc;
   %Save the estimation
   state_est(i,:) = ekf.state;
   vis.update(); 
end


%% Finally lets compare the original joint angle q and state prediction at
%the elbow
figure(2);clf;
plot(t,q);
title('Estimated and Actual joint angle');
hold on
%Recall our model has 3 prismatic joint for positioning, then 3 revolute
%shoulder joints and 1 elbow joint (the 7th joint)
plot(t,state_est(:,7));
legend('Actual','Estimated')
%Since our marker noise is super low estimation is basically pefect.


figure(3);clf;
plot(t,iteration_time*1000);
title(['Iteration Time, mean = ' num2str(mean(iteration_time)*1000) ' ms']);
xlabel('time (s)');
ylabel('Time per iteration (ms)');




