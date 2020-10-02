%Add main folder
addpath('..\');
%Create Model and visualize
mdl = rlCModel('..\Models\Lower_Body_Torso.xml');

mdl.forwardPosition();
mdl.base = 'rankle3';
%mdl.g = [0 0 0]';
mdl.forwardPosition();

sens_frames = {'rankle0', 'lankle0', 'rknee0', 'lknee0', 'mid_asis'};
%sens_frames = {'lankle0'};
for i=1:numel(sens_frames)
   sens = SensorCore([sens_frames{i} '_s']);
   sens.addDecorator('gyroscope');
   sens.addDecorator('accelerometer');
   T = eye(4);
   T(1:3,4) = [0 0 0]';
   mdl.addSensor(sens,sens_frames{i},T);
end

%Lets try to add another position sesnor at ankle
sens = SensorCore('pos_s');
sens.addDecorator('position');
mdl.addSensor(sens,'lankle3',eye(4));


mdl.forwardPosition();
%mdl.base = 'rankle3';
%mdl.forwardPosition();
vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();



dt = 0.01;
t = 0:dt:20;

q1 = cos(2*t);
dq1 = -2*sin(2*t);
ddq1 = -2*2*cos(2*t);

q2 = -q1;
dq2 = -dq1;
ddq2 = -ddq1;

%Run the trajectory and get measurements
mes = SensorMeasurement(numel(t),numel(mdl.sensors));

rknee_ind = ismember({mdl.joints.name},'rknee_j0');
lknee_ind = ismember({mdl.joints.name},'lknee_j0');


for i=1:numel(t)
    
    mdl.position(rknee_ind) = q1(i);
    mdl.velocity(rknee_ind) = dq1(i);
    mdl.acceleration(rknee_ind) = ddq1(i);
    mdl.position(lknee_ind) = q2(i);
    mdl.velocity(lknee_ind) = dq2(i);
    mdl.acceleration(lknee_ind) = ddq2(i);
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    
    mes(i,:) = SensorMeasurement(mdl.sensors);
    
   %vis.update();
end

%% Set up EKF
mdl.position(:) = 0;
mdl.velocity(:) = 0;
mdl.acceleration(:) = 0;
mdl.position(rknee_ind) = q1(1);
mdl.velocity(rknee_ind) = dq1(1);
mdl.acceleration(rknee_ind) = ddq1(1);
mdl.position(lknee_ind) = q2(1);
mdl.velocity(lknee_ind) = dq2(1);
mdl.acceleration(lknee_ind) = ddq2(1);
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

ekf = JC_EKF_Q_DQ_DDQ(mdl);

%Set up obs noise
ekf.observation_noise = diag(repmat([0.001 0.001 0.001 0.1 0.1 0.1],1,5));
%With additional position sensor
ekf.observation_noise = diag([repmat([0.001 0.001 0.001 0.1 0.1 0.1],1,5) 0.0001 0.0001 0.0001]);


G = [];
dim = numel(mdl.joints);
for i=ekf.sizeX/dim:-1:1
    G = [G; ones(dim,1)*dt^i / factorial(i)];
end
%Set up proc noise
P_tmp = G*G'*100;
P = zeros(size(P_tmp));
for i=1:ekf.sizeX/dim
    for j=i:ekf.sizeX/dim
        P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
            diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
    end
end

%What if there is no noise in position 
P = P+P' - diag(diag(P));

matches = [1:5]';
matches = [matches matches];
%% Do the estimation 

vel = zeros(numel(t),6);
accel = zeros(numel(t),3);
ekf_state = zeros(numel(t),ekf.sizeX);
mes_est = mes;

f_lankle = mdl.getFrameByName('lankle3');

g = mdl.g;

mes_i = SensorMeasurement(1,size(mes,2));
for i=1:size(mes,2)
mes_i(i).size = mes(1,i).size;
mes_i(i).type = mes(1,i).type;    
end

for i=1:numel(t)
    
    mes_i.setMesArray(mes(i,:).getMesArray());
    %mes_i(5).mes(3) = 0.1;
    %mes_i(2).mes(2) = 0.1;
    
    ekf.run_iteration(dt,mes_i,matches);
    ekf_state(i,:) = ekf.state;
    mes_est(i,:) = ekf.makeMeasure(ekf.state);
    vel(i,:) = f_lankle.v;
    accel(i,:) = f_lankle.t(1:3,1:3)*(f_lankle.a(4:6) + cross(f_lankle.v(1:3),f_lankle.v(4:6)))-g;
    vis.update();
end










