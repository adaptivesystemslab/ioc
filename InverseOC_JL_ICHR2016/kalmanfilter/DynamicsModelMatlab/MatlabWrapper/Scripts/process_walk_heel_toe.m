%% Walking using alternating constrants
addpath('../');

addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));
data_path = 'C:\aslab\data\walking\mocap\';
imu_dat_path = 'C:\aslab\data\walking\imu\';
mdlFilepath = '..\Models\avatar_v04_p_bothLegs_rev01.xml';
% 
% addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\common\ars'));
% data_path = 'D:\aslab\data\Walking_2019\Subject01\mocap\';
% imu_dat_path = 'D:\aslab\data\Walking_2019\Subject01\imu\';
% mdlFilepath = '..\Models\avatar_v04_p_bothLegs_rev01.xml';

visualize = 1;

%% Load the TRC File and build the model using Harrington 2007 for hip joint positions
time_stamp = '_2019_03_12_11_10_15';

%Load TRC data
trc = parseTRC([data_path 'templateafter1.trc']);
trcData = trc.data;
%Convert all to meters
m_names = fieldnames(trcData);
m_names = m_names(3:end);
for i=1:numel(m_names)
   trcData.(m_names{i}) = trcData.(m_names{i})/1000; 
end
%Mocap frames used for model building
m_indeces = 1:100;

%Load the model
mdl = rlCModel(mdlFilepath);

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
LL = (norm(rasis-rankle) + norm(lasis-lankle))/2;

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
mdl.position(4:6) = eul;
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
t_lankle_to_lheel = [lheel(1) lheel(2) 0]'-f_in.t(1:3,4);
t_lankle_to_ltoe = [ltoe(1) ltoe(2) 0]'-f_in.t(1:3,4);
t_lankle_to_lfoot_mid = [lfoot_mid(1) lfoot_mid(2) 0]'-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lankle_to_lheel;
%Right
t_indx = find(contains({mdl.transforms.name},'post:b_right_foot_to_pre:b_right_toe'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
rfoot_mid = (rtoe+rheel)/2;
t_rankle_to_rheel = [rheel(1) rheel(2) 0]'-f_in.t(1:3,4);
t_rankle_to_rtoe = [rtoe(1) rtoe(2) 0]'-f_in.t(1:3,4);
t_rankle_to_rfoot_mid = [rfoot_mid(1) rfoot_mid(2) 0]'-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_rankle_to_rheel;
mdl.forwardPosition();

%% Attach IMUS
imu_indeces = m_indeces;

imu_attachement = {...
    {'post:b_pelvis','TORSO_IMU_L','TORSO_IMU_T','TORSO_IMU_R'};...
    {'post:b_right_upperleg','KNEE_R_IMU_L','KNEE_R_IMU_T','KNEE_R_IMU_R'};...
    {'post:b_left_upperleg','KNEE_L_IMU_L','KNEE_L_IMU_T','KNEE_L_IMU_R'};...
    {'post:b_right_calf','ANKLE_R_IMU_L','ANKLE_R_IMU_T','ANKLE_R_IMU_R'};...
    {'post:b_left_calf','ANKLE_L_IMU_L','ANKLE_L_IMU_T','ANKLE_L_IMU_R'};...
    {'post:b_right_foot','FOOT_R_IMU_L','FOOT_R_IMU_T','FOOT_R_IMU_R'};...
    {'post:b_left_foot','FOOT_L_IMU_L','FOOT_L_IMU_T','FOOT_L_IMU_R'};...
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
f = mdl.getFrameByName('post:b_right_upperleg');
f_to_b = SE3.fastInverse(f.t)*b.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl.addSensor(sens,f.name,T);

sens = SensorCore('l_yaw');
sens.addDecorator('yaw');
sens.base = b;
f = mdl.getFrameByName('post:b_left_upperleg');
f_to_b = SE3.fastInverse(f.t)*b.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl.addSensor(sens,f.name,T);

sens = SensorCore('t_yaw');
sens.addDecorator('yaw');
f = mdl.getFrameByName('post:b_pelvis:rotateZ');
T = eye(4);
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl.addSensor(sens,f.name,T);

%% Add zero accelerationfor prismatic sensor
f = mdl.getFrameByName('post:b_pelvis:prismZ');
sens = SensorCore('prism_accel');
sens.addDecorator('accelerometer');
T = eye(4);
T(1:3,1:3) = f.t(1:3,1:3)';
T(1:3,4) = f.t(1:3,1:3)'*[0.1 0 0]';
mdl.addSensor(sens,'post:b_pelvis:prismZ',T);


%% Load IMU DATA
imu_files = {['0006667D7185' time_stamp '.header'],...
    ['00066681EB52' time_stamp '.header'],...
    ['0006667D713B' time_stamp '.header'],...
    ['00066681EA9B' time_stamp '.header'],...
    ['0006667D713A' time_stamp '.header'],...
    ['00066681EA98' time_stamp '.header'],...
    ['00066681EA9C' time_stamp '.header'],...
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
    
    accel = lowpass(imus(i).accelerometerCalibrated(1:min_samples,:),15,40);
    gyro = lowpass(imus(i).gyroscopeCalibrated(1:min_samples,:),15,40);
    
    mes = [mes gyro accel];
end

%Add yaw measurement
num_yaw = sum(contains({mdl.sensors.type},'yaw'));
mes = [mes zeros(size(mes,1),num_yaw)];

%Add prismatic acceleration to 0 measurement
mes = [mes repmat(mdl.g',size(mes,1),1)];

%Subtract off manually figured out X drift in TORSO sensor
%mes(:,1) = mes(:,1)-0.0036;

%% Set up and run EKF 
mdl.position(:) = 0;mdl.velocity(:) = 0;mdl.acceleration(:) = 0;
mdl.forwardPosition();mdl.forwardVelocity();mdl.forwardAcceleration();

%EKF with constrained left leg
ekf_lh = MC_EKF_Q_DQ_DDQ();
ekf_lh.dt = 1/imus(1).samplingRate;
ekf_lh.addModel(mdl);
load P_init.mat
ekf_lh.covariance = P_init;
ekf_lh.observation_noise = diag(repmat([0.01 0.01 0.01 1 1 1],1,7));
%Add yaw sensor noise
ekf_lh.observation_noise = blkdiag(ekf_lh.observation_noise,0.5*eye(num_yaw));
%Add acceleration to zero noise
ekf_lh.observation_noise = blkdiag(ekf_lh.observation_noise,0.1*eye(3));


%ekf_lc.process_noise(4:6,4:6) = ekf_lc.process_noise(4:6,4:6)+1e-6*eye(3);
%Add the left heel at current position constraint
fc_lh = mdl.getFrameByName('post:b_left_toe');
T = mdl.getFrameByName('post:b_left_toe').t;
C_lh = EKF_Constraint_BFConst('LC',mdl,'post:b_left_toe',eye(4),T);
C_lh.type(:) = false; C_lh.type([4 5 6]) = true;
C_lhv = EKF_Constraint_BFConst_vel('LCv',mdl,'post:b_left_toe',eye(4),zeros(6,1));
C_lhv.type(:) = false; C_lhv.type([4 5 6]) = true;
ekf_lh.addConstraint(C_lh);

%EKF with constrained right heel
ekf_rh = MC_EKF_Q_DQ_DDQ();
ekf_rh.dt = 1/imus(1).samplingRate;
ekf_rh.addModel(mdl);
ekf_rh.covariance = ekf_lh.covariance;
ekf_rh.observation_noise = ekf_lh.observation_noise;
ekf_rh.process_noise = ekf_lh.process_noise;
%Add the left foot at current position constraint
T = mdl.getFrameByName('post:b_right_toe').t;
fc_rh = mdl.getFrameByName('post:b_right_toe');
C_rh = EKF_Constraint_BFConst('RC',mdl,'post:b_right_toe',eye(4),T);
C_rh.type(:) = false; C_rh.type([4 5 6]) = true;
C_rhv = EKF_Constraint_BFConst_vel('RCv',mdl,'post:b_right_toe',eye(4),zeros(6,1));
C_rhv.type(:) = false; C_rv.type([4 5 6]) = true;
ekf_rh.addConstraint(C_rh);


ekf_rt = MC_EKF_Q_DQ_DDQ();
ekf_rt.dt = 1/imus(1).samplingRate;
ekf_rt.addModel(mdl);
ekf_rt.covariance = ekf_lh.covariance;
ekf_rt.observation_noise = ekf_lh.observation_noise;
ekf_rt.process_noise = ekf_lh.process_noise;
%Add right toe constraint
f = mdl.getFrameByName('post:b_right_toe');
%Place the constraint at the toe 
T= eye(4);
T(1:3,4) = f.t(1:3,1:3)'*([rtoe(1:2)-f.t(1:2,4); 0]);
C_rt = EKF_Constraint_BFConst('RT',mdl,'post:b_right_toe',T,eye(4));
C_rt.type(:) = false; C_rt.type([4 5 6]) = true;
ekf_rt.addConstraint(C_rt);

%Left Toe
ekf_lt = MC_EKF_Q_DQ_DDQ();
ekf_lt.dt = 1/imus(1).samplingRate;
ekf_lt.addModel(mdl);
ekf_lt.covariance = ekf_lh.covariance;
ekf_lt.observation_noise = ekf_lh.observation_noise;
ekf_lt.process_noise = ekf_lh.process_noise;
%Add right toe constraint
f = mdl.getFrameByName('post:b_left_toe');
%Place the constraint at the toe 
T= eye(4);
T(1:3,4) = f.t(1:3,1:3)'*([ltoe(1:2)-f.t(1:2,4); 0]);
C_lt = EKF_Constraint_BFConst('LT',mdl,'post:b_left_toe',T,eye(4));
C_lt.type(:) = false; C_lt.type([4 5 6]) = true;
ekf_lt.addConstraint(C_lt);



%% Visualize
if visualize
vis = rlVisualizer('vis',640,480);
mdl.forwardPosition();

vis.addModel(mdl);
vis.update();
end

%% Run EKF and see what happens

%Cur active left or right
cur_active = 'lh';
cur_count = 0;

ekf_state = zeros(size(mes,1),ekf_lh.sizeX);

eN_lh = zeros(size(mes,1),1);
eN_rh = zeros(size(mes,1),1);
eS_lh = zeros(size(mes,1),1);
eSN_lh = zeros(size(mes,1),1);
eP_lh = zeros(size(mes,1),1);
eS_rh = zeros(size(mes,1),1);
eSN_rh = zeros(size(mes,1),1);
eP_rh = zeros(size(mes,1),1);

eN_rt = zeros(size(mes,1),1);
eS_rt = zeros(size(mes,1),1);
eSN_rt = zeros(size(mes,1),1);
eP_rt = zeros(size(mes,1),1);

eN_lt = zeros(size(mes,1),1);
eS_lt = zeros(size(mes,1),1);
eSN_lt = zeros(size(mes,1),1);
eP_lt = zeros(size(mes,1),1);

%Selected foot -1: left 1:right
cur_foot = zeros(size(mes,1),1);
z_rh = zeros(size(mes,1),1);
z_lh = zeros(size(mes,1),1);
mes_pred = mes;
%Velocity of the feet
vel_l = zeros(size(mes,1),3);
vel_r = zeros(size(mes,1),3);
%Positions of the constrained feet
foot_pos = zeros(size(mes,1),2);

%Use yaw sensors
%num_yaw = 0;

%Number of samples that must have ghigher mean error
mean_samples = 1;

dt = 1/imus(1).samplingRate;
for i=3:size(mes,1)
    z = mes(i,:)';
    %Run Left Heel Constrained
    ekf_lh.run_iteration(dt,z);
    z_rh(i) = fc_rh.t(3,4);
    
    %Run Right Constrained
    ekf_rh.run_iteration(dt,z);
    z_lh(i) = fc_lh.t(3,4);
    %Run right toe constrained
    ekf_rt.run_iteration(dt,z);
    %Run left toe constrained
    ekf_lt.run_iteration(dt,z);
    
    %Re-run measurement prediction
    z_hat_lh = ekf_lh.makeMeasure(ekf_lh.state);
    %Get feet velocities from left constrained ekf
    if(strcmp(cur_active,'lh'))
        vel_l(i,:) = fc_lh.v(4:6);
        vel_r(i,:) = fc_rh.v(4:6);
    end
    H_lh = ekf_lh.makeH(ekf_lh.state); 
    S_lh = H_lh*ekf_lh.covariance*H_lh'+ekf_lh.observation_noise;
    S_lh = S_lh(1:end-num_yaw,1:end-num_yaw);
    S_lh_i = inv(S_lh);
    z_myaw = z(1:end-num_yaw);
    z_hat_lh_myaw = z_hat_lh(1:end-num_yaw);
    eS_lh(i) = (z_myaw-z_hat_lh_myaw)'*S_lh_i*(z_myaw-z_hat_lh_myaw);
    eSN_lh(i) = eS_lh(i)/norm(S_lh_i);
    eP_lh(i) = (z_myaw-z_hat_lh_myaw)'/(ekf_lh.inov_covariance(1:end-num_yaw,1:end-num_yaw))*(z_myaw-z_hat_lh_myaw);
    eN_lh(i) = (z_myaw-z_hat_lh_myaw)'*(z_myaw-z_hat_lh_myaw);
    
    z_hat_rh = ekf_rh.makeMeasure(ekf_rh.state);
    %Get feet velocities from left constrained ekf
    if(strcmp(cur_active,'rh'))
        vel_l(i,:) = fc_lh.v(4:6);
        vel_r(i,:) = fc_rh.v(4:6);
    end
    H_rh = ekf_rh.makeH(ekf_rh.state);
    S_rh = H_rh*ekf_rh.covariance*H_rh'+ekf_rh.observation_noise;
    S_rh = S_rh(1:end-num_yaw,1:end-num_yaw);
    S_rh_i = inv(S_rh);
    z_hat_rh_myaw = z_hat_rh(1:end-num_yaw);
    eS_rh(i) = (z_myaw-z_hat_rh_myaw)'*S_rh_i*(z_myaw-z_hat_rh_myaw);
    eSN_rh(i) = eS_rh(i)/norm(S_rh_i);
    eP_rh(i) = (z_myaw-z_hat_rh_myaw)'/(ekf_rh.inov_covariance(1:end-num_yaw,1:end-num_yaw))*(z_myaw-z_hat_rh_myaw);
    eN_rh(i) = (z_myaw-z_hat_rh_myaw)'*(z_myaw-z_hat_rh_myaw);
    
    z_hat_rt = ekf_rt.makeMeasure(ekf_rt.state);
    if(strcmp(cur_active,'rt'))
        vel_l(i,:) = fc_lh.v(4:6);
        vel_r(i,:) = fc_rh.v(4:6);
    end
    H_rt = ekf_rt.makeH(ekf_rt.state);
    S_rt = H_rt*ekf_rt.covariance*H_rt'+ekf_rt.observation_noise;
    S_rt = S_rt(1:end-num_yaw,1:end-num_yaw);
    S_rt_i = inv(S_rt);
    z_hat_rt_myaw = z_hat_rt(1:end-num_yaw);
    eS_rt(i) = (z_myaw-z_hat_rt_myaw)'*S_rt_i*(z_myaw-z_hat_rt_myaw);
    eSN_rt(i) = eS_rt(i)/norm(S_rt_i);
    eP_rt(i) = (z_myaw-z_hat_rt_myaw)'/(ekf_rt.inov_covariance(1:end-num_yaw,1:end-num_yaw))*(z_myaw-z_hat_rt_myaw);
    eN_rt(i) = (z_myaw-z_hat_rt_myaw)'*(z_myaw-z_hat_rt_myaw);  
    
    
    z_hat_lt = ekf_lt.makeMeasure(ekf_lt.state);
    if(strcmp(cur_active,'lt'))
        vel_l(i,:) = fc_lh.v(4:6);
        vel_r(i,:) = fc_rh.v(4:6);
    end
    H_lt = ekf_lt.makeH(ekf_lt.state);
    S_lt = H_lt*ekf_lt.covariance*H_lt'+ekf_lt.observation_noise;
    S_lt = S_lt(1:end-num_yaw,1:end-num_yaw);
    S_lt_i = inv(S_lt);
    z_hat_lt_myaw = z_hat_lt(1:end-num_yaw);
    eS_lt(i) = (z_myaw-z_hat_lt_myaw)'*S_lt_i*(z_myaw-z_hat_lt_myaw);
    eSN_lt(i) = eS_lt(i)/norm(S_lt_i);
    eP_lt(i) = (z_myaw-z_hat_lt_myaw)'/(ekf_lt.inov_covariance(1:end-num_yaw,1:end-num_yaw))*(z_myaw-z_hat_lt_myaw);
    eN_lt(i) = (z_myaw-z_hat_lt_myaw)'*(z_myaw-z_hat_lt_myaw);
    
    %Look at the error difference and make sure right leg isnt swinging.
    %if( (cur_active == 'l' && (mean(eSN_l(i-1:i)-eSN_r(i-1:i))) > 0.01 && z(6*4-4)>0.8 && z(6*5-4)>0.8 && cur_count > 15) || ...
    %        (cur_active == 'l' && cur_count > 15 && z_right(i) < -0.1))
    if(     (strcmp(cur_active,'lh') && (mean(eSN_lh(i-mean_samples:i)-eSN_lt(i-mean_samples:i))) > 0 && norm(vel_r(i,:)) < 1.5 && cur_count > 1))
        %We think left leg is swinging, but left leg error is larger than
        %right by a lot, switch constraint to right leg
        cur_active = 'lt';
        cur_count = 0;
    elseif( (strcmp(cur_active,'lh') && (mean(eSN_lh(i-mean_samples:i)-eSN_rh(i-mean_samples:i))) > 0 && norm(vel_r(i,:)) < 1.5 && cur_count > 10))
        cur_active = 'rh';
        cur_count = 0;
    elseif( (strcmp(cur_active,'lt') && (mean(eSN_lt(i-mean_samples:i)-eSN_rh(i-mean_samples:i))) > 0 && norm(vel_r(i,:)) < 1.5 && cur_count > 1))
        cur_active = 'rh';
        cur_count = 0;
    elseif( (strcmp(cur_active,'rh') && (mean(eSN_rh(i-mean_samples:i)-eSN_rt(i-mean_samples:i))) > 0 && norm(vel_l(i,:)) < 1.5 && cur_count > 1))
        cur_active = 'rt';
        cur_count = 0;
    elseif( (strcmp(cur_active,'rh') && (mean(eSN_rh(i-mean_samples:i)-eSN_lh(i-mean_samples:i))) > 0 && norm(vel_l(i,:)) < 1.5 && cur_count > 10))
        cur_active = 'lh';
        cur_count = 0;
    elseif( (strcmp(cur_active,'rt') && (mean(eSN_rt(i-mean_samples:i)-eSN_lh(i-mean_samples:i))) > 0 && norm(vel_l(i,:)) < 1.5 && cur_count > 1))
        cur_active = 'lh';
        cur_count = 0;
    end
    
    if(strcmp(cur_active,'lh'))
        cur_foot(i) = -1;
        mes_pred(i,:)=ekf_lh.makeMeasure(ekf_lh.state);
        ekf_rh.state = ekf_lh.state;
        ekf_rh.covariance = ekf_lh.covariance;
        ekf_rt.state = ekf_lh.state;
        ekf_rt.covariance = ekf_lh.covariance;
        ekf_lt.state = ekf_lh.state;
        ekf_lt.covariance = ekf_lh.covariance;
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        C_lt.Tw(1:3,4) = [C_lt.sensors(1).transform(1:2,4); 0];
        C_rh.Tw(1:3,4) = [C_rh.sensors(1).transform(1:2,4); 0];
        C_rt.Tw(1:3,4) = [C_rt.sensors(1).transform(1:2,4); 0];
        
        foot_pos(i,:) = fc_lh.t(1:2,4);
    elseif strcmp(cur_active,'rh')
        cur_foot(i) = 1;
        mes_pred(i,:)=ekf_rh.makeMeasure(ekf_rh.state);
        ekf_lh.state = ekf_rh.state;
        ekf_lh.covariance = ekf_rh.covariance;
        ekf_rt.state = ekf_rh.state;
        ekf_rt.covariance = ekf_rh.covariance;
        ekf_lt.state = ekf_rh.state;
        ekf_lt.covariance = ekf_rh.covariance;
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        C_lh.Tw(1:3,4) = [C_lh.sensors(1).transform(1:2,4); 0];
        C_lt.Tw(1:3,4) = [C_lt.sensors(1).transform(1:2,4); 0];
        C_rt.Tw(1:3,4) = [C_rt.sensors(1).transform(1:2,4); 0];
        foot_pos(i,:) = fc_rh.t(1:2,4);
    elseif strcmp(cur_active,'rt')
        cur_foot(i) = 2;
        mes_pred(i,:)=ekf_rt.makeMeasure(ekf_rt.state);
        ekf_lh.state = ekf_rt.state;
        ekf_lh.covariance = ekf_rt.covariance;
        ekf_rh.state = ekf_rt.state;
        ekf_rh.covariance = ekf_rt.covariance;
        ekf_lt.state = ekf_rt.state;
        ekf_lt.covariance = ekf_rt.covariance;
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        C_lh.Tw(1:3,4) = [C_lh.sensors(1).transform(1:2,4); 0];
        C_lt.Tw(1:3,4) = [C_lt.sensors(1).transform(1:2,4); 0];
        C_rh.Tw(1:3,4) = [C_rh.sensors(1).transform(1:2,4); 0];
        foot_pos(i,:) = fc_rh.t(1:2,4);
    elseif strcmp(cur_active,'lt')
        cur_foot(i) = -2;
        mes_pred(i,:)=ekf_lt.makeMeasure(ekf_lt.state);
        ekf_lh.state = ekf_lt.state;
        ekf_lh.covariance = ekf_lt.covariance;
        ekf_rh.state = ekf_lt.state;
        ekf_rh.covariance = ekf_lt.covariance;
        ekf_rt.state = ekf_lt.state;
        ekf_rt.covariance = ekf_lt.covariance;
        %Update constraints by fixing toes to Z = 0 at the current X and Y
        C_lh.Tw(1:3,4) = [C_lh.sensors(1).transform(1:2,4); 0];
        C_rt.Tw(1:3,4) = [C_rt.sensors(1).transform(1:2,4); 0];
        C_rh.Tw(1:3,4) = [C_rh.sensors(1).transform(1:2,4); 0];
        foot_pos(i,:) = fc_lh.t(1:2,4);
    end
    
    ekf_state(i,:) = ekf_lh.state;
    cur_count = cur_count+1;
    
    disp(['Frame: ' num2str(i) '/' num2str(size(mes,1)) ', Foot: ' cur_active ', el-er: ' num2str(eN_lh(i)-eN_rh(i))]);
    if mod(i,50) == 0 && visualize
        vis.update();
    end
    
    if mod(i,10) == 0 && i > 200 && visualize
        vis.addMarker(num2str(i),[foot_pos(i,:) 0]');
    end
end

figure(1);
plot(ekf_state(:, 2));


%% Look at specific peices

indx_start =2700;

for i=indx_start:size(mes,1)
   mdl.position = ekf_state(i,1:numel(mdl.joints));
   mdl.position(1:3) = 0;
   mdl.forwardPosition();
   vis.update();
   pause(dt);
end



