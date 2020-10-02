%% Walking using alternating constrants
addpath('../');

addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));
data_path = 'C:\aslab\data\walking\Subject03\mocap\';
imu_dat_path = 'C:\aslab\data\walking\Subject03\imu\';
mdlFilepath = '..\Models\avatar_v04_p_bothLegs_zyx_hip_no_ank.xml';
% 
% addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\common\ars'));
% data_path = 'D:\aslab\data\Walking_2019\Subject01\mocap\';
% imu_dat_path = 'D:\aslab\data\Walking_2019\Subject01\imu\';
% mdlFilepath = '..\Models\avatar_v04_p_bothLegs_rev01.xml';

visualize = 1;

%% Load the TRC File and build the model using Harrington 2007 for hip joint positions
time_stamp = '_2019_03_28_14_59_58';

%Load TRC data
trc = parseTRC([data_path 'templatebefore1_smoothed.trc']);
trcData = trc.data;
%Convert all to meters
m_names = fieldnames(trcData);
m_names = m_names(3:end);
for i=1:numel(m_names)
   trcData.(m_names{i}) = trcData.(m_names{i})/1000; 
end
%Mocap frames used for model building
m_indeces = 1000:1500;

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

%% Attach IMUS
imu_indeces = m_indeces;

imu_attachement = {...
    {'post:b_pelvis','TORSO_IMU_L','TORSO_IMU_T','TORSO_IMU_R'};...
    {'post:b_right_upperleg','KNEE_R_IMU_L','KNEE_R_IMU_T','KNEE_R_IMU_R'};...
    {'post:b_left_upperleg','KNEE_L_IMU_L','KNEE_L_IMU_T','KNEE_L_IMU_R'};...
    {'post:b_right_calf','ANKLE_R_IMU_L','ANKLE_R_IMU_T','ANKLE_R_IMU_R'};...
    {'post:b_left_calf','ANKLE_L_IMU_L','ANKLE_L_IMU_T','ANKLE_L_IMU_R'};...
    ...%{'post:b_right_foot','FOOT_R_IMU_L','FOOT_R_IMU_T','FOOT_R_IMU_R'};...
    ...%{'post:b_left_foot','FOOT_L_IMU_L','FOOT_L_IMU_T','FOOT_L_IMU_R'};...
    };

for i=1:numel(imu_attachement)
    f_name =  imu_attachement{i}{1};
    frame = mdl.getFrameByName(f_name);
    
    P = mean(trcData.(imu_attachement{i}{2})(imu_indeces,:));
    Q = mean(trcData.(imu_attachement{i}{3})(imu_indeces,:));
    R = mean(trcData.(imu_attachement{i}{4})(imu_indeces,:));
    [R_imu_0, R_0_imu] = points2rot(P,Q,R);
    t_0 =  mean([P;Q;R])' - frame.t(1:3,4);
    t_f = frame.t(1:3,1:3)'*t_0;
    
    sens = SensorCore(imu_attachement{i}{2}(1:end-2));
    sens.addDecorator('gyroscope');
    sens.addDecorator('accelerometer');
    T = eye(4);
    T(1:3,1:3) = frame.t(1:3,1:3)'*R_0_imu';
    T(1:3,4) = t_f;
    mdl.addSensor(sens,frame.name,T)
end
mdl.forwardPosition();

%% Attach YAW sensors to each leg 
b = mdl.getFrameByName('post:b_pelvis:rotateZ');

sens = SensorCore('r_yaw');
sens.addDecorator('yaw');
sens.base = b;
f = mdl.getFrameByName('post:b_right_upperleg:rotateZ');
f_to_b = SE3.fastInverse(f.t)*b.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl.addSensor(sens,f.name,T);

sens = SensorCore('l_yaw');
sens.addDecorator('yaw');
sens.base = b;
f = mdl.getFrameByName('post:b_left_upperleg:rotateZ');
f_to_b = SE3.fastInverse(f.t)*b.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl.addSensor(sens,f.name,T);

% sens = SensorCore('t_yaw');
% sens.addDecorator('yaw');
% f = mdl.getFrameByName('post:b_pelvis:rotateZ');
% T = eye(4);
% T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
% mdl.addSensor(sens,f.name,T);

%% Add zero accelerationfor prismatic sensor
f = mdl.getFrameByName('post:b_pelvis:prismZ');
sens = SensorCore('prism_accel');
sens.addDecorator('accelerometer');
T = eye(4);
T(1:3,1:3) = f.t(1:3,1:3)';
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
%mdl.addSensor(sens,'post:b_pelvis:prismZ',T);


%% Load IMU DATA
imu_files = {['0006667D7185' time_stamp '.header'],...
    ['0006667D713B' time_stamp '.header'],...
    ['00066681EA9C' time_stamp '.header'],...
    ['0006667D713A' time_stamp '.header'],...
    ['00066681EA9B' time_stamp '.header'],...
    ...%['0006667BA9F2' time_stamp '.header'],...
    ...%['00066681EA98' time_stamp '.header'],...
    };
imus = tinySensorDataHandle.empty();
num_samples = [];
for i=1:numel(imu_files)
    imu = arsLoader([imu_dat_path imu_files{i}]);
    imus(i) = imu;
    num_samples(i) = size(imu.accelerometerCalibrated,1);
end

min_samples = min([num_samples]);

mes = [];
for i = 1:numel(imus)
    
    %accel = lowpass(imus(i).accelerometerCalibrated(1:min_samples,:),100,imu(1).samplingRate);
    %gyro = lowpass(imus(i).gyroscopeCalibrated(1:min_samples,:),100,imu(1).samplingRate);
    
    gyro_int = cumtrapz(imus(i).time-imus(i).time(1),imus(i).gyroscopeCalibrated);
    gyro_bias = (gyro_int(end,:)-gyro_int(1,:))/(imus(i).time(end)-imus(i).time(1));
    gyro_unbias = imus(i).gyroscopeCalibrated-...
        repmat(gyro_bias,size(imus(i).gyroscopeCalibrated,1),1);
    
    
    gyro= imus(i).gyroscopeCalibrated(1:min_samples,:);
    accel =imus(i).accelerometerCalibrated(1:min_samples,:);
    mes = [mes gyro accel];
end

%Add yaw measurement
num_yaw = sum(contains({mdl.sensors.type},'yaw'));
num_imu = sum(contains({mdl.sensors.type},'gyroscope'));
mes = [mes zeros(size(mes,1),num_yaw)];

%Add prismatic acceleration to 0 measurement
%mes = [mes repmat(mdl.g',size(mes,1),1)];


%% Set up and run EKF 
mdl.position(:) = 0;mdl.velocity(:) = 0;mdl.acceleration(:) = 0;
mdl.forwardPosition();mdl.forwardVelocity();mdl.forwardAcceleration();

%mdl.transforms(1).t(1:2,4) = [10 , -10]';
%mdl.position(1:2) = [-10, 10]';
%mdl.forwardPosition();

%EKF with constrained left leg
ekf_lc = MC_EKF_Q_DQ_DDQ();
ekf_lc.dt = 1/imus(1).samplingRate;
ekf_lc.addModel(mdl);
%load P_init.mat
ekf_lc.covariance = ekf_lc.process_noise;
ekf_lc.observation_noise = diag(repmat([0.01 0.01 0.01 0.1 0.1 0.1],1,num_imu));
%Add yaw sensor noise
ekf_lc.observation_noise = blkdiag(ekf_lc.observation_noise,0.01*eye(num_yaw));
%Add acceleration to zero noise
%ekf_lc.observation_noise = blkdiag(ekf_lc.observation_noise,1000*eye(3));


%ekf_lc.process_noise(4:6,4:6) = ekf_lc.process_noise(4:6,4:6)+1e-6*eye(3);
%Add the left foot at current position constraint
fc_l = mdl.getFrameByName('post:b_left_toe');
T = mdl.getFrameByName('post:b_left_toe').t;
C_l = EKF_Constraint_BFConst('LC',mdl,'post:b_left_toe',eye(4),T);
C_l.type(:) = false; C_l.type([4 5 6]) = true;
C_lv = EKF_Constraint_BFConst_vel('LCv',mdl,'post:b_left_toe',eye(4),zeros(6,1));
C_lv.type(:) = false; C_lv.type([4 5 6]) = true;
ekf_lc.addConstraint(C_l);
%ekf_lc.addConstraint(C_lv);

%EKF with constrained right leg
ekf_rc = MC_EKF_Q_DQ_DDQ();
ekf_rc.dt = 1/imus(1).samplingRate;
ekf_rc.addModel(mdl);
ekf_rc.covariance = ekf_lc.covariance;
ekf_rc.observation_noise = ekf_lc.observation_noise;
ekf_rc.process_noise = ekf_lc.process_noise;
%Add the left foot at current position constraint
T = mdl.getFrameByName('post:b_right_toe').t;
fc_r = mdl.getFrameByName('post:b_right_toe');
C_r = EKF_Constraint_BFConst('RC',mdl,'post:b_right_toe',eye(4),T);
C_r.type(:) = false; C_r.type([4 5 6]) = true;
C_rv = EKF_Constraint_BFConst_vel('RCv',mdl,'post:b_right_toe',eye(4),zeros(6,1));
C_rv.type(:) = false; C_rv.type([4 5 6]) = true;
ekf_rc.addConstraint(C_r);
%ekf_rc.addConstraint(C_rv);

%% Visualize
if visualize
vis = rlVisualizer('vis',640,480);
mdl.forwardPosition();

vis.addModel(mdl);
vis.update();
end

%% Run EKF and see what happens
ekf_state = zeros(size(mes,1),ekf_lc.sizeX);
%Selected foot -1: left 1:right
cur_foot = zeros(size(mes,1),1);
z_right = zeros(size(mes,1),1);
z_left = zeros(size(mes,1),1);
mes_pred = mes;
%Velocity of the feet
vel_l = zeros(size(mes,1),3);
vel_r = zeros(size(mes,1),3);
%Positions of the constrained feet
foot_pos = zeros(size(mes,1),2);

%This is the look ahead stuff
looking = false;    %This flag becomes true after we load up the first horizon
horizon = 4;       %How far we look ahead
%Saved States
hor_state_l = zeros(horizon,ekf_lc.sizeX);
hor_state_r = zeros(horizon,ekf_rc.sizeX);
%Saved Covariances
hor_cov_l = zeros(ekf_lc.sizeX,ekf_lc.sizeX,horizon);
hor_cov_r = zeros(ekf_rc.sizeX,ekf_rc.sizeX,horizon);
%Saved projected covariances
hor_icov_l = zeros(ekf_lc.sizeZ,ekf_lc.sizeZ,horizon);
hor_icov_r = zeros(ekf_rc.sizeZ,ekf_lc.sizeZ,horizon);
%Jacobians
hor_H_l = zeros(ekf_lc.sizeZ,ekf_lc.sizeX,horizon);
hor_H_r = zeros(ekf_rc.sizeZ,ekf_rc.sizeX,horizon);
%Saved predicted measurements
hor_mes_l = zeros(horizon,ekf_lc.sizeZ);
hor_mes_r = zeros(horizon,ekf_rc.sizeZ);
%The errors
hor_Es_l = zeros(horizon,1);
hor_Es_r = zeros(horizon,1);

%Preload the look ahead
start_indx = 680;
for i=start_indx:start_indx+horizon-1
    indx = i-start_indx+1;
    z = mes(i,:)';
    z_myaw = z(1:end-num_yaw);
    
    %Run Left Constrained and populate horizon vars
    ekf_lc.run_iteration(dt,z);
    vel_r(i,:) = fc_r.v(4:6);
    hor_state_l(indx,:) = ekf_lc.state;
    hor_cov_l(:,:,indx) = ekf_lc.covariance;
    %Re-run mes prediction
    hor_mes_l(indx,:) = ekf_lc.makeMeasure(ekf_lc.state);
    z_hat_left_myaw = hor_mes_l(indx,1:end-num_yaw)';
    H_left = ekf_lc.makeH(ekf_lc.state); 
    hor_H_l(:,:,indx) = H_left;
    hor_icov_l(:,:,indx) = H_left*ekf_lc.covariance*H_left'+ekf_lc.observation_noise;
    S_left = hor_icov_l(1:end-num_yaw,1:end-num_yaw,indx);
    S_left_i = inv(S_left);
    hor_Es_l(indx) = (z_myaw-z_hat_left_myaw)'*S_left_i*(z_myaw-z_hat_left_myaw);
    
    %Run Right Constrained and populate horizon vars
    ekf_rc.run_iteration(dt,z);
    vel_l(i,:) = fc_l.v(4:6);
    hor_state_r(indx,:) = ekf_rc.state;
    hor_cov_r(:,:,indx) = ekf_rc.covariance;
    %Re-run mes prediction
    hor_mes_r(indx,:) = ekf_rc.makeMeasure(ekf_rc.state);
    z_hat_right_myaw = hor_mes_r(indx,1:end-num_yaw)';
    H_right = ekf_rc.makeH(ekf_rc.state); 
    hor_H_r(:,:,indx) = H_right;
    hor_icov_r(:,:,indx) = H_right*ekf_rc.covariance*H_right'+ekf_rc.observation_noise;
    S_right = hor_icov_r(1:end-num_yaw,1:end-num_yaw,indx);
    S_right_i = inv(S_right);
    hor_Es_r(indx) = (z_myaw-z_hat_right_myaw)'*S_right_i*(z_myaw-z_hat_right_myaw);
end
%% Decide on the starting foot
El = sum(hor_Es_l)/horizon;
Er = sum(hor_Es_r)/horizon;

if El < Er
   cur_active = 'l'; 
else
   cur_active = 'r';
end

%% Run the rest with looking ahead
for i=start_indx+1:size(mes,1)-horizon
    
    %Pull out the actual measurement horizon
    z_hor = mes(i:i+horizon-1,:);
    
    %For the active foot we just need to run the next iteration 
    if cur_active == 'r'
        %Save the current estiamted state with right foot on ground
        cur_state = hor_state_r(1,:)';
        cur_cov = hor_cov_r(:,:,1);
        z_myaw = z_hor(end,1:end-num_yaw)';
        %Run Right Constrained for the last measurement
        ekf_rc.state = hor_state_r(end,:)';
        ekf_rc.covariance = hor_cov_r(:,:,end);
        ekf_rc.run_iteration(dt,z_hor(end,:)');
        vel_l(i+horizon-1,:) = fc_l.v(4:6);
        vel_r(i+horizon-1,:) = fc_r.v(4:6);
        %Shift all the things down and put in the last one
        hor_state_r(1:end-1) = hor_state_r(2:end);
        hor_state_r(end,:) = ekf_rc.state;
        hor_cov_r(:,:,1:end-1) = hor_cov_r(:,:,2:end);
        hor_cov_r(:,:,end) = ekf_rc.covariance;
        %Re-run mes prediction
        hor_mes_r(1:end-1,:) = hor_mes_r(2:end,:);
        hor_mes_r(end,:) = ekf_rc.makeMeasure(ekf_rc.state);
        z_hat_right_myaw = hor_mes_r(end,1:end-num_yaw)';
        H_right = ekf_rc.makeH(ekf_rc.state); 
        hor_H_r(:,:,1:end-1) = hor_H_r(:,:,2:end);
        hor_H_r(:,:,end) = H_right;
        hor_icov_r(:,:,1:end-1) = hor_icov_r(:,:,2:end);
        hor_icov_r(:,:,end) = H_right*ekf_rc.covariance*H_right'+ekf_rc.observation_noise;
        S_right = hor_icov_r(1:end-num_yaw,1:end-num_yaw,end);
        S_right_i = inv(S_right);
        hor_Es_r(1:end-1) = hor_Es_r(2:end);
        hor_Es_r(end) = (z_myaw-z_hat_right_myaw)'*S_right_i*(z_myaw-z_hat_right_myaw);
        
        %Now for the left foot we have to run the entire horizon starting
        %at the first right constrained state
        ekf_lc.state = cur_state;
        ekf_lc.covariance = cur_cov;
        %Run for the entire horizon assuming left foot is constrainted
        for indx = 1:size(z_hor,1)
            %This is the measurement we are looking at
            z = z_hor(indx,:)';
            z_myaw = z(1:end-num_yaw);
            %Run Left Constrained and populate horizon vars
            ekf_lc.run_iteration(dt,z);
            hor_state_l(indx,:) = ekf_lc.state;
            hor_cov_l(:,:,indx) = ekf_lc.covariance;
            %Re-run mes prediction
            hor_mes_l(indx,:) = ekf_lc.makeMeasure(ekf_lc.state);
            z_hat_left_myaw = hor_mes_l(indx,1:end-num_yaw)';
            H_left = ekf_lc.makeH(ekf_lc.state); 
            hor_H_l(:,:,indx) = H_left;
            hor_icov_l(:,:,indx) = H_left*ekf_lc.covariance*H_left'+ekf_lc.observation_noise;
            S_left = hor_icov_l(1:end-num_yaw,1:end-num_yaw,indx);
            S_left_i = inv(S_left);
            hor_Es_l(indx) = (z_myaw-z_hat_left_myaw)'*S_left_i*(z_myaw-z_hat_left_myaw);
        end
    else
        %Save the current estiamted state with right foot on ground
        cur_state = hor_state_l(1,:)';
        cur_cov = hor_cov_l(:,:,1);
        z_myaw = z_hor(end,1:end-num_yaw)';
        %Run Left Constrained for the last measurement
        ekf_lc.state = hor_state_l(end,:)';
        ekf_lc.covariance = hor_cov_l(:,:,end);
        ekf_lc.run_iteration(dt,z_hor(end,:)');
        vel_l(i+horizon-1,:) = fc_l.v(4:6);
        vel_r(i+horizon-1,:) = fc_r.v(4:6);
        %Shift all the things down and put in the last one
        hor_state_l(1:end-1) = hor_state_l(2:end);
        hor_state_l(end,:) = ekf_lc.state;
        hor_cov_l(:,:,1:end-1) = hor_cov_l(:,:,2:end);
        hor_cov_l(:,:,end) = ekf_lc.covariance;
        %Re-run mes prediction
        hor_mes_l(1:end-1,:) = hor_mes_l(2:end,:);
        hor_mes_l(end,:) = ekf_lc.makeMeasure(ekf_lc.state);
        z_hat_left_myaw = hor_mes_l(end,1:end-num_yaw)';
        H_left = ekf_lc.makeH(ekf_lc.state); 
        hor_H_l(:,:,1:end-1) = hor_H_l(:,:,2:end);
        hor_H_l(:,:,end) = H_left;
        hor_icov_l(:,:,1:end-1) = hor_icov_l(:,:,2:end);
        hor_icov_l(:,:,end) = H_left*ekf_lc.covariance*H_left'+ekf_lc.observation_noise;
        S_left = hor_icov_r(1:end-num_yaw,1:end-num_yaw,end);
        S_left_i = inv(S_left);
        hor_Es_l(1:end-1) = hor_Es_l(2:end);
        hor_Es_l(end) = (z_myaw-z_hat_left_myaw)'*S_left_i*(z_myaw-z_hat_left_myaw);
        
        %Now for the left foot we have to run the entire horizon starting
        %at the first right constrained state
        ekf_rc.state = cur_state;
        ekf_rc.covariance = cur_cov;
        %Run for the entire horizon assuming left foot is constrainted
        for indx = 1:size(z_hor,1)
            %This is the measurement we are looking at
            z = z_hor(indx,:)';
            z_myaw = z(1:end-num_yaw);
            %Run Left Constrained and populate horizon vars
            ekf_rc.run_iteration(dt,z);
            hor_state_r(indx,:) = ekf_rc.state;
            hor_cov_r(:,:,indx) = ekf_rc.covariance;
            %Re-run mes prediction
            hor_mes_r(indx,:) = ekf_rc.makeMeasure(ekf_rc.state);
            z_hat_right_myaw = hor_mes_r(indx,1:end-num_yaw)';
            H_right = ekf_rc.makeH(ekf_rc.state); 
            hor_H_r(:,:,indx) = H_right;
            hor_icov_r(:,:,indx) = H_right*ekf_rc.covariance*H_right'+ekf_rc.observation_noise;
            S_right = hor_icov_r(1:end-num_yaw,1:end-num_yaw,indx);
            S_right_i = inv(S_right);
            hor_Es_r(indx) = (z_myaw-z_hat_right_myaw)'*S_right_i*(z_myaw-z_hat_right_myaw);
        end       
    end
    
    %Look at the error difference over the horizon
    Er = sum(hor_Es_r)/horizon;
    El = sum(hor_Es_l)/horizon;
    is_r = isempty(find(hor_Es_l - hor_Es_r < 0,1));
    is_l = isempty(find(hor_Es_r - hor_Es_l < 0,1));
    if cur_active == 'l' && Er < El && norm(vel_r(i,:)) < 1.4
        %We think left leg is swinging, but left leg error is larger than
        %right by a lot, switch constraint to right leg
        cur_active = 'r';
    elseif cur_active == 'r' && El < Er && norm(vel_l(i,:)) < 1.4
        cur_active = 'l';
        cur_count = 0;
    end
    
    if(cur_active == 'l')
        cur_foot(i) = -1;
        ekf_state(i,:) = hor_state_l(1,:);
        mes_pred(i,:)=ekf_lc.makeMeasure(hor_state_l(1,:)');
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        T = fc_r.t;
        T(3,4) = 0;
        C_r.Tw = T;
        foot_pos(i,:) = fc_l.t(1:2,4);
    else
        cur_foot(i) = 1;
        ekf_state(i,:) = hor_state_r(1,:);
        mes_pred(i,:)=ekf_lc.makeMeasure(hor_state_r(1,:)');
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        T = fc_l.t;
        T(3,4) = 0;
        C_l.Tw = T;
        foot_pos(i,:) = fc_r.t(1:2,4);
    end
    
    disp(['Frame: ' num2str(i) '/' num2str(size(mes,1)) ', Foot: ' cur_active ', el-er: ' num2str(Er-El)]);
    if visualize
        vis.update();
    end
    if mod(i,30) == 0 && i > 200 && visualize
        vis.addMarker(num2str(i),[foot_pos(i,:) 0]');
    end
end


%% Pull out step length

%Pull out foot positions when right foot is on the ground
foot_pos_r = foot_pos(cur_foot==1,:);
%Pull out foot positions when left foot is on the ground
foot_pos_l = foot_pos(cur_foot==-1,:);

%foot sections 
sections_r = abs(diff(foot_pos_r(:,1)))+abs(diff(foot_pos_r(:,2))) > 0.1;
sections_r_indxs = find(sections_r);
sections_l = abs(diff(foot_pos_l(:,1)))+abs(diff(foot_pos_l(:,2))) > 0.1;
sections_l_indxs = find(sections_l);

%Compute mean foot positions and step lengths
mean_foot_pos_r =zeros(numel(sections_r_indxs),2);
step_lengths_r = zeros(numel(sections_r_indxs),1);

for i=1:numel(sections_r_indxs)-1
    section = foot_pos_r(sections_r_indxs(i):sections_r_indxs(i+1),:);
    mean_foot_pos_r(i,:) = mean(section);
end

for i=1:size(mean_foot_pos_r,1)-2
   step_lengths_r(i) = norm(mean_foot_pos_r(i+1,:)-mean_foot_pos_r(i,:));
end

mean_foot_pos_l =zeros(numel(sections_l_indxs),2);
step_lengths_l = zeros(numel(sections_l_indxs),1);

for i=1:numel(sections_l_indxs)-1
    section = foot_pos_l(sections_l_indxs(i):sections_l_indxs(i+1),:);
    mean_foot_pos_l(i,:) = mean(section);
end

for i=1:size(mean_foot_pos_l,1)-2
   step_lengths_l(i) = norm(mean_foot_pos_l(i+1,:)-mean_foot_pos_l(i,:));
end

