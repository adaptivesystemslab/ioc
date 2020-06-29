%% Walking using alternating constrants
addpath('../');

addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));
data_path = 'C:\aslab\data\walking\RobohubWalk\Subject01\';
mdlFilepath = '..\Models\avatar_v04_p_bothLegs_zyx_hip.xml';

mdlTRCPath = [data_path 'static_before03.trc'];
ikTRCPath = [data_path 'circle__CW_sr01.trc'];

visualize = 1;
%% Load the TRC File and build the model using Harrington 2007 for hip joint positions

%Load TRC data
trc = parseTRC(mdlTRCPath);
trcData = trc.data;
%Convert all to meters
m_names = fieldnames(trcData);
m_names = m_names(3:end);
for i=1:numel(m_names)
   trcData.(m_names{i}) = (rotx(pi/2)*trcData.(m_names{i})'/1000)'; 
end
%Mocap frames used for model building
m_indeces = 100:300;

%Load the model
mdl = rlCModel(mdlFilepath);
mdl.forwardPosition
vis = rlVisualizer('vis',640,480)
vis.addModel(mdl)
vis.update


%Rotation calculate based on the mean of mid->front and right for y
%Here we get the world to base transform to match Mocap Data
front = mean((trcData.ASIS_R(m_indeces,:)+trcData.ASIS_L(m_indeces,:))/2);
back = mean((trcData.PSIS_R(m_indeces,:)+trcData.PSIS_L(m_indeces,:))/2);
mid = (front+back)/2;
left = mean(trcData.ASIS_L(m_indeces,:));
right = mean(trcData.ASIS_R(m_indeces,:));

%This is rotation from world to body orieantation in mocap
[~,R] = points2rot([front(1:2) 0],[mid(1:2) 0],[left(1:2) 0]);

lasis = mean(R*trcData.ASIS_L(m_indeces,:)',2);
rasis = mean(R*trcData.ASIS_R(m_indeces,:)',2);
lback = mean(R*trcData.PSIS_L(m_indeces,:)',2);
rback = mean(R*trcData.PSIS_R(m_indeces,:)',2);
lknee = (mean(R*trcData.KNEE_L_LAT(m_indeces,:)',2) + mean(R*trcData.KNEE_L_MED(m_indeces,:)',2))/2;
rknee = (mean(R*trcData.KNEE_R_LAT(m_indeces,:)',2) + mean(R*trcData.KNEE_R_MED(m_indeces,:)',2))/2;
lankle = (mean(R*trcData.ANKLE_L_LAT(m_indeces,:)',2) + mean(R*trcData.ANKLE_L_MED(m_indeces,:)',2))/2;
rankle = (mean(R*trcData.ANKLE_R_LAT(m_indeces,:)',2) + mean(R*trcData.ANKLE_R_MED(m_indeces,:)',2))/2;
ltoe = (mean(R*trcData.FOOT_L_LAT(m_indeces,:)',2) + mean(R*trcData.FOOT_L_MED(m_indeces,:)',2))/2;
rtoe = (mean(R*trcData.FOOT_R_LAT(m_indeces,:)',2) + mean(R*trcData.FOOT_R_MED(m_indeces,:)',2))/2;
lheel = mean(R*trcData.HEEL_L(m_indeces,:)',2);
rheel = mean(R*trcData.HEEL_R(m_indeces,:)',2);

%Distance between hips as calculated in

%@ARTICLE{harrington2007prediction,
%    author = {Harrington, ME and Zavatsky, AB and Lawson, SEM and Yuan, Z and Theologis,TN},
%    title = {Prediction of the hip joint centre in adults, children, and patients
%    with cerebral palsy based on magnetic resonance imaging},
%    journal = {Journal of biomechanics},
%    year = {2007},
%    volume = {40}
%}

%Pelvis Width calculated using x and y only
PW = norm(lasis(1:2) - rasis(1:2));
%Pelvis Depth calculated using x and y only
PD = norm(abs((lasis(1:2)+rasis(1:2))/2 - (lback(1:2)+rback(1:2))/2));
%Leg Length
LL = (norm(rasis-rheel) + norm(lasis-lheel))/2;

%So going from the middle the joint center is predicted in mm as
x = -0.24*PD-9.9/1000;
y = 0.28*PD+0.16*PW+7.9/1000;
z = -0.16*PW-0.04*LL-7.1/1000;
% Where x is from middle to front, y is from middle to side, and z is
% up
asismid = (lasis + rasis)/2;
hipcentre = (lasis + rasis + lback + rback) / 4;
hipcentre_to_asismid = asismid - hipcentre;

%Set model's origin to middle of pelvis
mdl.transforms(1).t = eye(4);
mdl.transforms(1).t(1:3,4) = hipcentre;

%Set model's mid pelvis orientation
eul = rotm2eul(R','XYZ');
mdl.position(4:6) = -eul;
mdl.forwardPosition();

%Set the transform from middle of pelvis to hip joint centers
%Note that hip joint centers start at post:b_pelvis body while the
%harrington x y z is in x = front, y = side, z = up frames so we rotate it
%since out model is y forward x side z up

%Here we have pelvis middle to hip joint centers in world frame
%hipcentre to left hip joint center
hipcentre_to_lhipjc = R*([x y z]' + hipcentre_to_asismid);
%hipcentre to right hip joint center
hipcentre_to_rhipjc = R*([x -y z]' + hipcentre_to_asismid);

%We convert it to our model frame 
%Set the transforms
t_indx = find(contains({mdl.transforms.name},'post:b_pelvis_to_pre:b_left_upperleg'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
trans.t(1:3,4) = f_in.t(1:3,1:3)'*hipcentre_to_lhipjc;

t_indx = find(contains({mdl.transforms.name},'post:b_pelvis_to_pre:b_right_upperleg'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
trans.t(1:3,4) = f_in.t(1:3,1:3)'*hipcentre_to_rhipjc;
mdl.forwardPosition();

%Hip to knee
%Left
t_indx = find(contains({mdl.transforms.name},'post:b_left_upperleg_to_pre:b_left_calf'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
t_lhip_to_lknee = lknee-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lhip_to_lknee;
%Right
t_indx = find(contains({mdl.transforms.name},'post:b_right_upperleg_to_pre:b_right_calf'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
t_rhip_to_rknee = rknee-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_rhip_to_rknee;
mdl.forwardPosition();

%Knee to Ankle
%Left
t_indx = find(contains({mdl.transforms.name},'post:b_left_calf_to_pre:b_left_foot'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
t_lknee_to_lankle = lankle-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lknee_to_lankle;
%Right
t_indx = find(contains({mdl.transforms.name},'post:b_right_calf_to_pre:b_right_foot'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
t_lknee_to_lankle = rankle-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lknee_to_lankle;
mdl.forwardPosition();

%Ankle to middle of foot
%Left
t_indx = find(contains({mdl.transforms.name},'post:b_left_foot_to_pre:b_left_toe'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
lfoot_mid = (ltoe+lheel)/2;
t_lankle_to_ltoe = [lheel(1) lheel(2) 0]'-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lankle_to_ltoe;
%Right
t_indx = find(contains({mdl.transforms.name},'post:b_right_foot_to_pre:b_right_toe'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
rfoot_mid = (rtoe+rheel)/2;
t_rankle_to_rtoe = [rheel(1) rheel(2) 0]'-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_rankle_to_rtoe;
mdl.forwardPosition();

%% Attach Markers
marker_attachement = {...
    {'post:b_pelvis','ASIS_R','ASIS_L','PSIS_R','PSIS_L'};...
    {'post:b_right_upperleg','KNEE_R_MED','KNEE_R_LAT'};...
    {'post:b_left_upperleg','KNEE_L_MED','KNEE_L_LAT'};...
    {'post:b_right_calf','ANKLE_R_LAT','ANKLE_R_MED'};...
    {'post:b_left_calf','ANKLE_L_LAT','ANKLE_L_MED'};...
    {'post:b_right_foot','FOOT_R_LAT','FOOT_R_MED','HEEL_R'};...
    {'post:b_left_foot','FOOT_L_LAT','FOOT_L_MED','HEEL_L'};...
    };

for i=1:numel(marker_attachement)
    f_name =  marker_attachement{i}{1};
    frame = mdl.getFrameByName(f_name);
    
    for j=2:numel(marker_attachement{i})
        m_name = marker_attachement{i}{j};
        
        
        %Figure out the mean position
        marker_pos = trcData.(m_name);
        %Rotate to body frame
        mean_m_pos = R*mean(marker_pos(m_indeces,:))';
        
        %Create the sensor
        sens = SensorCore(marker_attachement{i}{j});
        sens.addDecorator('position');
        
        %Attach the sensor to the model
        t_w = mean_m_pos-frame.t(1:3,4); %Translation in world frame
        t_f = frame.t(1:3,1:3)'*t_w;     %Translation in frame
        T = eye(4);
        T(1:3,4) = t_f;
        mdl.addSensor(sens,frame.name,T)
    end
end
mdl.forwardPosition();

%% Visualize the model
if visualize
   vis = rlVisualizer('vis',640,480);
   vis.addModel(mdl);
   vis.update();
end

%% Load the TRC that will be used for IK
ik_trc = parseTRC(ikTRCPath);
ik_trcData = ik_trc.data;
%Convert all to meters
m_names = fieldnames(ik_trcData);
m_names = m_names(3:end);
for i=1:numel(m_names)
   ik_trcData.(m_names{i}) = (rotx(pi/2)*ik_trcData.(m_names{i})'/1000)'; 
end

%% Position the model such that it matches the first IK TRC frame
first_indeces = 1:5;
front = mean((ik_trcData.ASIS_R(first_indeces,:)+ik_trcData.ASIS_L(first_indeces,:))/2);
back = mean((ik_trcData.PSIS_R(first_indeces,:)+ik_trcData.PSIS_L(first_indeces,:))/2);
mid = (front+back)/2;
left = mean(ik_trcData.ASIS_L(first_indeces,:));
right = mean(ik_trcData.ASIS_R(first_indeces,:));
%This is rotation from world to body orieantation in mocap
[~,R] = points2rot([front(1:2) 0],[mid(1:2) 0],[left(1:2) 0]);
asismid = (lasis + rasis)/2;
hipcentre = mid;
%Set model's origin to middle of pelvis
mdl.transforms(1).t = eye(4);
mdl.transforms(1).t(1:3,4) = hipcentre;
%Set model's mid pelvis orientation
eul = rotm2eul(R','XYZ');
mdl.position(4:6) = -eul;
mdl.forwardPosition();

%Add actual markers to see if we got it close or not
for i=1:numel(mdl.sensors)
    m_name = mdl.sensors(i).name;
    vis.addMarker(m_name,mean(ik_trcData.(m_name)(first_indeces,:)),[1 1 0 1]); 
end
vis.update();

%% Build measurement vectors
mes = [];
for i=1:numel(marker_attachement)
    for j = 2:numel(marker_attachement{i})
        mes = [mes ik_trcData.(marker_attachement{i}{j})];
    end
end
mes_obj = SensorMeasurement(size(mes,1),numel(mdl.sensors));
[mes_obj.type] = deal(1);
[mes_obj.size] = deal(3);
mes_obj.setMesArray(mes);

%% Set up magik EKF for IK 

ekf = EKF_Q_DQ_DDQ(mdl);
ekf.observation_noise = ekf.observation_noise*0.01;
eta = 10;
dt = 1/trc.DataRate;
%Process noise
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
ekf_rarm.process_noise = P;

%% Run EKF to converge to first frame

matches = [1:numel(mdl.sensors)]';
matches = [matches matches];
for i=1:100
    z = mes_obj(1,:);    
    ekf.run_iteration(dt,z,matches);
end
vis.update();

%% Run ekf for all frames
state_est = zeros(size(mes,1),numel(ekf.state));
for i=1:size(mes,1)
    
    %This is the measurement for current timestep, must be vertical vector
    z = mes_obj(i,:);    
    ekf.run_iteration(dt,z);
    state_est(i,:) = ekf.state;
    %Draw markers
    %for j=1:numel(m_names)
    %    vis.addMarker(m_names{j},trc.data.(m_names{j})(i,:));
    %end
    disp(['Frame: ' num2str(i) ' / ' num2str(size(mes,1))]);
    if mod(i,10) == 0
        vis.update
    end
end

%% Low pass Q 
lpFilt = designfilt('lowpassfir','PassbandFrequency',12, ...
         'StopbandFrequency',14,'PassbandRipple',0.5, ...
         'StopbandAttenuation',65,'DesignMethod','kaiserwin','SampleRate',trc.DataRate);

state_est_filt = state_est;
dof = numel(mdl.joints);
for i=1:numel(mdl.joints)
   state_est_filt(:,i) = filtfilt(lpFilt,state_est(:,i)); 
   state_est_filt(:,i+dof) = [diff(state_est_filt(:,i))/dt ; 0];
   state_est_filt(:,i+dof*2) = [diff(state_est_filt(:,i+dof))/dt; 0];
end

%% Create a new model with IMUs attached to it
mdl_imu = rlCModel(mdlFilepath);
mdl_imu.forwardPosition();

%Make this new model same as trc model
trc_transforms = {mdl.transforms.name};
for i=1:numel(mdl_imu.transforms)
    if(strcmp(mdl_imu.transforms(i).type,'fixed'))
        t_indx = find(contains(trc_transforms,mdl_imu.transforms(i).name),1,'first');
        mdl_imu.transforms(i).t = mdl.transforms(t_indx).t;
    end
end
%Set to world origin 
mdl_imu.transforms(1).t = eye(4);
mdl_imu.forwardPosition();

%% Attach fake IMUs to the model 
imu_attachement = {...
    {'post:b_pelvis','TORSO_IMU',[-0.1,0,0.1]'};...
    {'post:b_right_upperleg','KNEE_R',[0.05,0,-0.2]'};...
    {'post:b_left_upperleg','KNEE_L',[0.05,0,-0.2]'};...
    {'post:b_right_calf','ANKLE_R',[0.05,0,-0.2]'};...
    {'post:b_left_calf','ANKLE_L',[0.05,0,-0.2]'};...
    {'post:b_right_foot','FOOT_R',[0.1,0,0]'};...
    {'post:b_left_foot','FOOT_L',[0.1,0,0]'};...
    };
T  =eye(4);
for i=1:numel(imu_attachement)
   frame = mdl_imu.getFrameByName(imu_attachement{i}{1});
   sens = SensorCore(imu_attachement{i}{2});
   sens.addDecorator('accelerometer');
   sens.addDecorator('gyroscope');
   t_w = imu_attachement{i}{3}
   t_f = frame.t(1:3,1:3)'*t_w;
   T(1:3,1:3) = frame.t(1:3,1:3)';
   T(1:3,4) = t_f;
   mdl_imu.addSensor(sens,frame.name,T);
end
mdl_imu.forwardPosition();

%% Attach YAW sensors to each leg 
b = mdl_imu.getFrameByName('post:b_pelvis:rotateZ');

sens = SensorCore('r_yaw');
sens.addDecorator('yaw');
sens.base = b;
f = mdl_imu.getFrameByName('post:b_right_upperleg:rotateZ');
f_to_b = SE3.fastInverse(f.t)*b.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl_imu.addSensor(sens,f.name,T);

sens = SensorCore('l_yaw');
sens.addDecorator('yaw');
sens.base = b;
f = mdl_imu.getFrameByName('post:b_left_upperleg:rotateZ');
f_to_b = SE3.fastInverse(f.t)*b.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl_imu.addSensor(sens,f.name,T);

%% Place model at initial position
%Set model's origin to middle of pelvis
mdl_imu.transforms(1).t = eye(4);
mdl_imu.transforms(1).t(1:3,4) = hipcentre;

%Set model's mid pelvis orientation
mdl_imu.position(4:6) = -eul;
mdl_imu.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl_imu);
vis.update();

%% Generate IMU measurements
dof = numel(mdl_imu.joints);
z = vertcat(mdl_imu.sensors.measurement);
mes_imu = zeros(size(mes,1),numel(z));
for i=1:size(mes_imu,1)
   mdl_imu.position = state_est_filt(i,1:dof); 
   mdl_imu.velocity = state_est_filt(i,dof+1:dof*2); 
   mdl_imu.acceleration = state_est_filt(i,dof*2+1:end); 
   mdl_imu.forwardPosition;
   mdl_imu.forwardVelocity;
   mdl_imu.forwardAcceleration;
   
   mes_imu(i,:) = vertcat(mdl_imu.sensors.measurement);
   if mod(i,10)==0
   vis.update();
   end
end

%% Add noise to IMU measurements
ACCEL_NOISE = 0.1;
GYRO_NOSIE = 0.01;
YAW_NOISE = 0.1;

num_imu = sum(contains({mdl_imu.sensors.type},'accelerometer'));
num_yaw = sum(contains({mdl_imu.sensors.type},'yaw'));

mes_imu_noisy = mes_imu;
for i=1:num_imu
    mes_imu_noisy(:,i*6-5:i*6) = mes_imu_noisy(:,i*6-5:i*6)...
        + [ACCEL_NOISE*randn(size(mes,1),3) GYRO_NOSIE*randn(size(mes,1),3)];
end

%Make zero YAW
mes_imu_noisy(:,end-num_yaw+1:end) = 0;

%% Now reset the model and set up the alternating constrained EKF
mdl_imu.position = state_est_filt(1,1:dof);
mdl_imu.velocity = state_est_filt(1,dof+1:dof*2);
mdl_imu.acceleration = state_est_filt(1,dof*2+1:end);
mdl_imu.forwardPosition;
mdl_imu.forwardVelocity;
mdl_imu.forwardAcceleration;

%EKF with constrained left leg
ekf_lc = MC_EKF_Q_DQ_DDQ();
ekf_lc.dt = dt;
ekf_lc.addModel(mdl_imu);
%load P_init.mat
ekf_lc.covariance = ekf_lc.process_noise;
                                        %ACCEL        GYRO
ekf_lc.observation_noise = diag(repmat([0.1 0.1 0.1 0.01 0.01 0.01],1,7));
%Yaw noise
num_yaw = sum(contains({mdl_imu.sensors.type},'yaw'));
ekf_lc.observation_noise = blkdiag(ekf_lc.observation_noise,0.001*eye(num_yaw));


%Add the left foot at current position constraint
fc_l = mdl_imu.getFrameByName('post:b_left_toe');
T = mdl_imu.getFrameByName('post:b_left_toe').t;
C_l = EKF_Constraint_BFConst('LC',mdl_imu,'post:b_left_toe',eye(4),T);
C_l.type(:) = false; C_l.type([4 5 6]) = true;
C_lv = EKF_Constraint_BFConst_vel('LCv',mdl_imu,'post:b_left_toe',eye(4),zeros(6,1));
C_lv.type(:) = false; C_lv.type([4 5 6]) = true;
ekf_lc.addConstraint(C_l);

%EKF with constrained right leg
ekf_rc = MC_EKF_Q_DQ_DDQ();
ekf_rc.dt = dt;
ekf_rc.addModel(mdl_imu);
ekf_rc.covariance = ekf_lc.covariance;
ekf_rc.observation_noise = ekf_lc.observation_noise;
ekf_rc.process_noise = ekf_lc.process_noise;
%Add the left foot at current position constraint
T = mdl_imu.getFrameByName('post:b_right_toe').t;
fc_r = mdl_imu.getFrameByName('post:b_right_toe');
C_r = EKF_Constraint_BFConst('RC',mdl_imu,'post:b_right_toe',eye(4),T);
C_r.type(:) = false; C_r.type([4 5 6]) = true;
C_rv = EKF_Constraint_BFConst_vel('RCv',mdl_imu,'post:b_right_toe',eye(4),zeros(6,1));
C_rv.type(:) = false; C_rv.type([4 5 6]) = true;
ekf_rc.addConstraint(C_r);

%% Visualize
vis = rlVisualizer('vis',640,480);
mdl_imu.forwardPosition();
vis.addModel(mdl_imu);
vis.update();

%% Run EKF and see what happens

%Cur active left or right
cur_active = 'l';
cur_count = 0;

ekf_state = zeros(size(mes,1),ekf_lc.sizeX);

eN_l = zeros(size(mes,1),1);
eN_r = zeros(size(mes,1),1);
eS_l = zeros(size(mes,1),1);
eSN_l = zeros(size(mes,1),1);
eP_l = zeros(size(mes,1),1);
eS_r = zeros(size(mes,1),1);
eSN_r = zeros(size(mes,1),1);
eP_r = zeros(size(mes,1),1);
%Selected foot -1: left 1:right
cur_foot = zeros(size(mes,1),1);
z_right = zeros(size(mes,1),1);
z_left = zeros(size(mes,1),1);
mes_pred = mes_imu;
%Velocity of the feet
vel_l = zeros(size(mes,1),3);
vel_r = zeros(size(mes,1),3);
%Positions of the constrained feet
foot_pos = zeros(size(mes,1),2);

for i=100:size(mes,1)
    z = mes_imu(i,:)';
    %Run Left Constrained
    ekf_lc.run_iteration(dt,z);
    z_right(i) = fc_r.t(3,4);
    %Run Right Constrained
    ekf_rc.run_iteration(dt,z);
    z_left(i) = fc_l.t(3,4);
    
    %Re-run measurement prediction
    z_hat_left = ekf_lc.makeMeasure(ekf_lc.state);
    %Get feet velocities from left constrained ekf
    if(cur_active == 'l')
        vel_l(i,:) = fc_l.v(4:6);
        vel_r(i,:) = fc_r.v(4:6);
    end
    
    H_left = ekf_lc.makeH(ekf_lc.state); 
    S_left = H_left*ekf_lc.covariance*H_left'+ekf_lc.observation_noise;
    S_left = S_left(1:end-num_yaw,1:end-num_yaw);
    S_left_i = inv(S_left);
    z_myaw = z(1:end-num_yaw);
    z_hat_left_myaw = z_hat_left(1:end-num_yaw);
    eS_l(i) = (z_myaw-z_hat_left_myaw)'*S_left_i*(z_myaw-z_hat_left_myaw);
    eSN_l(i) = eS_l(i)/norm(S_left_i);
    eP_l(i) = (z_myaw-z_hat_left_myaw)'/(ekf_lc.inov_covariance(1:end-num_yaw,1:end-num_yaw))*(z_myaw-z_hat_left_myaw);
    eN_l(i) = (z_myaw-z_hat_left_myaw)'*(z_myaw-z_hat_left_myaw);
    
    
    z_hat_right = ekf_rc.makeMeasure(ekf_rc.state);
    %Get feet velocities from left constrained ekf
    if(cur_active == 'r')
        vel_l(i,:) = fc_l.v(4:6);
        vel_r(i,:) = fc_r.v(4:6);
    end
    H_right = ekf_rc.makeH(ekf_rc.state);
    S_right = H_right*ekf_rc.covariance*H_right'+ekf_rc.observation_noise;
    S_right = S_right(1:end-num_yaw,1:end-num_yaw);
    S_right_i = inv(S_right);
    z_hat_right_myaw = z_hat_right(1:end-num_yaw);
    eS_r(i) = (z_myaw-z_hat_right_myaw)'*S_right_i*(z_myaw-z_hat_right_myaw);
    eSN_r(i) = eS_r(i)/norm(S_right_i);
    eP_r(i) = (z_myaw-z_hat_right_myaw)'/(ekf_rc.inov_covariance(1:end-num_yaw,1:end-num_yaw))*(z_myaw-z_hat_right_myaw);
    eN_r(i) = (z_myaw-z_hat_right_myaw)'*(z_myaw-z_hat_right_myaw);
    
    e_d = norm(eS_r(i)-eS_l(i));
    
    
    %Look at the error difference and make sure right leg isnt swinging.
    %if( (cur_active == 'l' && (mean(eSN_l(i-1:i)-eSN_r(i-1:i))) > 0.01 && z(6*4-4)>0.8 && z(6*5-4)>0.8 && cur_count > 15) || ...
    %        (cur_active == 'l' && cur_count > 15 && z_right(i) < -0.1))
    if((cur_active == 'l' && eSN_l(i)-eSN_r(i) > -0.1 && norm(vel_r(i,:)) < 0.5  && cur_count > 30))
        %We think left leg is swinging, but left leg error is larger than
        %right by a lot, switch constraint to right leg
        cur_active = 'r';
        cur_count = 0;
    %elseif( (cur_active == 'r' && (mean(eSN_r(i-1:i)-eSN_l(i-1:i))) > 0.01 && z(6*4-4)>0.8 && z(6*5-4)>0.8 && cur_count > 15) || ...
    %        (cur_active == 'r' && cur_count > 15 && z_left(i) < -0.1))
    elseif( (cur_active == 'r' && eSN_r(i)-eSN_l(i) > -0.1 && norm(vel_l(i,:)) < 0.5 && cur_count > 30))
        cur_active = 'l';
        cur_count = 0;
    end
    
    if(cur_active == 'l')
        cur_foot(i) = -1;
        mes_pred(i,:)=ekf_lc.makeMeasure(ekf_lc.state);
        ekf_rc.state = ekf_lc.state;
        ekf_rc.covariance = ekf_lc.covariance;
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        T = fc_r.t;
        T(3,4) = 0;
        C_r.Tw = T;
        foot_pos(i,:) = fc_l.t(1:2,4);
    else
        cur_foot(i) = 1;
        mes_pred(i,:)=ekf_rc.makeMeasure(ekf_rc.state);
        ekf_lc.state = ekf_rc.state;
        ekf_lc.covariance = ekf_rc.covariance;
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        T = fc_l.t;
        T(3,4) = 0;
        C_l.Tw = T;
        foot_pos(i,:) = fc_r.t(1:2,4);
    end
    
    ekf_state(i,:) = ekf_lc.state;
    cur_count = cur_count+1;
    
    disp(['Frame: ' num2str(i) '/' num2str(size(mes,1)) ', Foot: ' cur_active ', el-er: ' num2str(eN_l(i)-eN_r(i))]);
    if visualize
        vis.update();
    end
    if mod(i,30) == 0 && i > 200 && visualize
        vis.addMarker(num2str(i),[foot_pos(i,:) 0]');
    end
end


