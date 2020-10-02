%Test Measurement Constrained EKF

%Add main folder
addpath('..\');

%Create Model and visualize
mdl = rlCModel('simple_3dof_model.xml');

%add sensor to end effector
sens = SensorCore('sens1');
sens.addDecorator('position');
%sens.addDecorator('velocity');
mdl.addSensor(sens,'frame_ee',eye(4));
sens = SensorCore('sens2');
sens.addDecorator('position');
%sens.addDecorator('velocity');
mdl.addSensor(sens,'frame2',eye(4));
sens = SensorCore('sens3');
sens.addDecorator('position');
%sens.addDecorator('velocity');
mdl.addSensor(sens,'frame4',eye(4));

mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

%Generate Motion 
t = [0:0.01:pi]';
dt  =t(2);
q1 = t;
q2 = zeros(size(t));
q3 = zeros(size(t));

q = [q1 q2 q2];
dq = ones(size(q));
dq(:,[2 3]) = 0;
mes = SensorMeasurement(numel(t),numel(mdl.sensors));

for i=1:numel(t);
    mdl.position = q(i,:);
    mdl.velocity = dq(i,:);
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mes(i,:) = SensorMeasurement(mdl.sensors);

    %vis.update();
    %pause(t(2));
end

%% CONSTRAINED ESTIMATION
%We set up a wall at x = 0.5 and try to use SVD constrained EKF 

%Reset Model
mdl.position(:) = 0;

%Create EKF
ekf = SVD_EKF_Q_DQ(mdl);
ekf.observation_noise = diag(repmat([0.01 0.01 0.01],1,numel(mdl.sensors)));
%ekf.observation_noise = diag([0.01 0.01 0.01 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6]);
%Velocity noise
dim = numel(mdl.joints);
eta = 10;
G = [ones(dim,1)*dt^2/2*eta; ones(dim,1)*dt*eta];
P_tmp = G*G';
P = zeros(size(P_tmp));
for i=1:2
    for j=i:2
        P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
            diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
    end
end
P = P+P' - diag(diag(P));
ekf.process_noise = P;


%ekf.upper_mes_bound(1) = 0.5;

ekf.upper_state_bound(1) = pi/4;
%ekf.upper_state_bound(4:6) = 2;
%ekf.lower_state_bound(:) = -pi/4;

ekf_state_est = zeros(numel(t),numel(ekf.state));
ekf_mes_est = zeros(numel(t),numel(ekf.makeMeasure(ekf.state).getMesArray));
%Estimate
for i=1:numel(t);
    
    if i==206
       disp('bp'); 
    end
    
    ekf.run_iteration(dt,mes(i,:),[1 1;2 2;3 3]);
    ekf_state_est(i,:) = ekf.state;
    ekf_mes_est(i,:) = ekf.makeMeasure(ekf.state).getMesArray;
    
    for j=1:numel(mes(i,:))
        vis.addMarker(num2str(j),mes(i,j).mes);
    end
    vis.update();
    pause(dt);
end




