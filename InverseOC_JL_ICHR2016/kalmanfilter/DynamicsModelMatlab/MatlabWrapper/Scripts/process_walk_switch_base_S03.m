%% Walking using alternating constrants
addpath('../');

addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));
data_path = 'C:\aslab\data\walking\Subject03\mocap\';
imu_dat_path = 'C:\aslab\data\walking\Subject03\imu\';
mdlFilepath = '..\Models\avatar_v04_p_bothLegs_zyx_hip.xml';
%
% addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\common\ars'));
% data_path = 'D:\aslab\data\Walking_2019\Subject01\mocap\';
% imu_dat_path = 'D:\aslab\data\Walking_2019\Subject01\imu\';
% mdlFilepath = '..\Models\avatar_v04_p_bothLegs_rev01.xml';

visualize = 1;

%% Load the TRC File and build the model using Harrington 2007 for hip joint positions
time_stamp = '_2019_03_28_14_53_09';

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
t_lknee_to_lankle = [lheel(1) lheel(2) 0]'-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lknee_to_lankle;
%Right
t_indx = find(contains({mdl.transforms.name},'post:b_right_calf_to_pre:b_right_foot'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
t_lknee_to_lankle = [rheel(1) rheel(2) 0]'-f_in.t(1:3,4);
trans.t(1:3,4) = f_in.t(1:3,1:3)'*t_lknee_to_lankle;
mdl.forwardPosition();

%Ankle to middle of foot
%Left
t_indx = find(contains({mdl.transforms.name},'post:b_left_foot_to_pre:b_left_toe'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
trans.t(1:3,4) = 0;
%Right
t_indx = find(contains({mdl.transforms.name},'post:b_right_foot_to_pre:b_right_toe'),1);
trans = mdl.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;
trans.t(1:3,4) = 0;
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

%% Attach YAW sensors to right leg wrt pelvis
f_yaw_base = mdl.getFrameByName('pre:b_left_upperleg:rotateZ');
sens_yaw = SensorCore('s_yaw');
sens_yaw.addDecorator('yaw');
sens_yaw.base = f_yaw_base;
f_yaw_r = mdl.getFrameByName('post:b_right_upperleg:rotateZ');
f_to_b = SE3.fastInverse(f_yaw_r.t)*f_yaw_base.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f_yaw_r.t(1:3,1:3)'*[0.1 0 0]';
T_yaw_r = T;
mdl.addSensor(sens_yaw,f_yaw_r.name,T);

f_yaw_l = mdl.getFrameByName('post:b_left_upperleg:rotateZ');
f_to_b = SE3.fastInverse(f_yaw_l.t)*f_yaw_base.t;
T = eye(4);
T(1:3,1:3)=f_to_b(1:3,1:3);
T(1:3,4) = f_yaw_l.t(1:3,1:3)'*[0.1 0 0]';
T_yaw_l = T;

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
    
    %accel = lowpass(imus(i).accelerometerCalibrated(1:min_samples,:),5,imu(1).samplingRate);
    %gyro = lowpass(imus(i).gyroscopeCalibrated(1:min_samples,:),5,imu(1).samplingRate);
    accel = imus(i).accelerometerCalibrated(1:min_samples,:);
    gyro = imus(i).gyroscopeCalibrated(1:min_samples,:);
    imus(i).time = imus(i).time(1:min_samples);
    imus(i).gyroscopeCalibrated = gyro;
    imus(i).accelerometerCalibrated = accel;
    
    mes = [mes gyro accel];
end

%Add yaw measurement
num_yaw = sum(contains({mdl.sensors.type},'yaw'));
mes = [mes zeros(size(mes,1),num_yaw)];

%Subtract off manually figured out X drift in TORSO sensor
%mes(:,1) = mes(:,1)-0.0036;

%% Set up and run EKF
mdl.position(:) = 0;mdl.velocity(:) = 0;mdl.acceleration(:) = 0;
mdl.forwardPosition();mdl.forwardVelocity();mdl.forwardAcceleration();

%Set left foot as base for the model
mdl.base = 'post:b_left_toe';
%Create EKF with left foot support
ekf = MC_EKF_Q_DQ_DDQ();
ekf.dt = 1/imus(1).samplingRate;
ekf.addModel(mdl);
ekf.covariance = ekf.process_noise;
%Noise for torso sensor
ekf.observation_noise = diag([0.01 0.01 0.01 1 1 1]);
ekf.observation_noise = blkdiag(ekf.observation_noise,diag(repmat([0.01 0.01 0.01 1 1 1],1,4)));
%Add yaw sensor noise
ekf.observation_noise = blkdiag(ekf.observation_noise,0.001*eye(num_yaw));
%Create Right Foot projection onto ground constraint
C_r = EKF_Constraint_BFConst('RC',mdl,'post:b_right_toe',eye(4),T);
C_r.type(:) = false; C_r.type([4 5 6]) = true;
%Save foot frame
fc_r = mdl.getFrameByName('post:b_right_toe');
%Keep track of state_l and cov_l
state_l = ekf.state;
cov_l = ekf.covariance;
%Create the left foot projection onto ground constraint
%C_l = EKF_Constraint_BFConst('LC',mdl,'post:b_left_toe',eye(4),T);
C_l = EKF_Constraint_BFConst('LC',mdl,'post:b_left_toe',eye(4),T);
C_l.type(:) = false; C_l.type([4 5 6]) = true;
%Save foot frame
fc_l = mdl.getFrameByName('post:b_left_toe');
%Keep track of state_r and cov_r
state_r = ekf.state;
cov_r = ekf.covariance;

%% Visualize
if visualize
    vis = rlVisualizer('vis',640,480);
    mdl.forwardPosition();
    vis.addModel(mdl);
    vis.update();
end

%% Run EKF and see what happens

%Lets start standing on right foot, the model base is already right
cur_active = 'l';
cur_count = 0;

ekf_state = zeros(size(mes,1),ekf.sizeX);

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
mes_pred = mes;
%Velocity of the feet
vel_l = zeros(size(mes,1),3);
vel_r = zeros(size(mes,1),3);
%Positions of the constrained feet
foot_pos = zeros(size(mes,1),2);

%Set base to left foot
mdl.base = 'post:b_left_toe';

%Lankle joints 
lankle_joints = find(contains({mdl.joints.name},'pre:b_left_foot:rotateX_to_post:b_left_foot:rotateX') |...
    contains({mdl.joints.name},'pre:b_left_foot:rotateZ_to_post:b_left_foot:rotateZ'));
%rankle joints
rankle_joints = find(contains({mdl.joints.name},'pre:b_right_foot:rotateX_to_post:b_right_foot:rotateX') |...
    contains({mdl.joints.name},'pre:b_right_foot:rotateZ_to_post:b_right_foot:rotateZ'));


dt = 1/imus(1).samplingRate;
dof = numel(mdl.position);
for i=660:size(mes,1)
    z = mes(i,:)';
    
    %Save current state
    cur_state = ekf.state;
    
    if cur_active == 'l'
        cur_foot(i) = -1;
        %Save state
        ekf_state(i,:) = state_l;
        
        %Set right ankle to zero and not moving
        indxs_r = [rankle_joints dof+rankle_joints dof*2+rankle_joints];
        indxs_l = [lankle_joints dof+lankle_joints dof*2+lankle_joints];
        state_l(indxs_r) = 0;
        state_l(dof+rankle_joints) = state_l(dof+lankle_joints);
        state_l(dof*2+rankle_joints) = state_l(dof*2+lankle_joints);
        cov_l(indxs_r,:) = cov_l(indxs_l,:);
        cov_l(:,indxs_r) = cov_l(:,indxs_l);
        
        %Project right foot to ground and collapse the covariance
        T = fc_r.t;
        T(3,4) = 0;
        C_r.Tw = T;
        mdl.calculateSensorJacobians();
        Ac = C_r.A;
        b = C_r.b;
        A = zeros(size(Ac,1),ekf.sizeX);
        A(:,1:size(Ac,2)) = Ac;
        W_k_inv  = cov_l;
        Y = W_k_inv*A'/(A*W_k_inv*A' + eye(size(A,1))*1e-5);
        %Project State
        state_r_proj = state_l - Y*(A*state_l - b);
        %Update covariance
        cov_r_proj = (eye(size(cov_l)) - Y*A)*cov_l;
        
        %Run regular left EKF
        ekf.state = state_l;
        ekf.covariance = cov_l;
        ekf.run_iteration(dt,z);
        state_l = ekf.state;
        cov_l = ekf.covariance;
        z_hat_left = ekf.makeMeasure(ekf.state);
        mes_pred(i,:) = z_hat_left;
        H_left = ekf.makeH(ekf.state);
        vel_r(i,:) = fc_r.v(4:6);
        foot_pos(i,:) = fc_l.t(1:2,4);
        
        %Switch base and run projected left
        mdl.position = state_r_proj(1:dof);
        mdl.velocity = state_r_proj(dof+1:dof*2);
        mdl.acceleration = state_r_proj(dof*2+1:end);
        mdl.forwardPosition;
        mdl.base = 'post:b_right_toe';
        tmp_state_r_base = [mdl.position; mdl.velocity; mdl.acceleration];
        tmp_state_r_base(rankle_joints+dof) = state_l(lankle_joints+dof);
        tmp_state_r_base(rankle_joints+dof*2) = state_l(lankle_joints+dof*2);
        ekf.state = tmp_state_r_base;
        ekf.covariance = cov_r_proj;
        ekf.run_iteration(dt,z);
        state_r = ekf.state;
        cov_r = ekf.covariance;
        z_hat_right = ekf.makeMeasure(ekf.state);
        H_right = ekf.makeH(ekf.state);
        
        %Switch back to right
        mdl.position = tmp_state_r_base(1:dof);
        mdl.velocity = tmp_state_r_base(dof+1:dof*2);
        mdl.acceleration = tmp_state_r_base(dof*2+1:end);
        mdl.forwardPosition;
        mdl.base = 'post:b_left_toe';
        %Update to latest right base results
        mdl.position = state_l(1:dof);
        mdl.velocity = state_l(dof+1:dof*2);
        mdl.acceleration = state_l(dof*2+1:end);
        mdl.forwardPosition;mdl.forwardVelocity;mdl.forwardAcceleration;
        
    elseif cur_active == 'r'
        cur_foot(i) = 1;
        %Save state
        ekf_state(i,:) = -state_r;
        
        %Set right ankle to zero and not moving
        indxs_r = [rankle_joints dof+rankle_joints dof*2+rankle_joints];
        indxs_l = [lankle_joints dof+lankle_joints dof*2+lankle_joints];
        state_r(indxs_l) = 0;
        state_r(dof+lankle_joints) = state_r(dof+rankle_joints);
        state_r(dof*2+lankle_joints) = state_r(dof*2+rankle_joints);
        cov_r(indxs_l,:) = cov_r(indxs_r,:);
        cov_r(:,indxs_l) = cov_r(:,indxs_r);
        
        %Project left foot to ground and collapse the covariance
        T = fc_l.t;
        T(3,4) = 0;
        C_l.Tw = T;
        mdl.calculateSensorJacobians();
        Ac = C_l.A;
        b = C_l.b;
        A = zeros(size(Ac,1),ekf.sizeX);
        A(:,1:size(Ac,2)) = Ac;
        W_k_inv  = cov_r;
        Y = W_k_inv*A'/(A*W_k_inv*A' + eye(size(A,1))*1e-5);
        %Project State
        state_l_proj = state_r - Y*(A*state_r - b);
        %Update covariance
        cov_l_proj = (eye(size(cov_r)) - Y*A)*cov_r;
        
        %Run regular right EKF
        ekf.state = state_r;
        ekf.covariance = cov_r;
        ekf.run_iteration(dt,z);
        state_r = ekf.state;
        cov_r = ekf.covariance;
        z_hat_right = ekf.makeMeasure(ekf.state);
        mes_pred(i,:) = z_hat_right;
        H_right = ekf.makeH(ekf.state);
        vel_l(i,:) = fc_l.v(4:6);
        foot_pos(i,:) = fc_r.t(1:2,4);
        
        %Switch base and run projected left
        mdl.position = state_l_proj(1:dof);
        mdl.velocity = state_l_proj(dof+1:dof*2);
        mdl.acceleration = state_l_proj(dof*2+1:end);
        mdl.forwardPosition;
        mdl.base = 'post:b_left_toe';
        tmp_state_l_base = [mdl.position; mdl.velocity; mdl.acceleration];
        tmp_state_l_base(lankle_joints+dof) = state_r(rankle_joints+dof);
        tmp_state_l_base(lankle_joints+dof*2) = state_r(rankle_joints+dof*2);
        ekf.state = tmp_state_l_base;
        ekf.covariance = cov_l_proj;
        ekf.run_iteration(dt,z);
        state_l = ekf.state;
        cov_l = ekf.covariance;
        z_hat_left = ekf.makeMeasure(ekf.state);
        H_left = ekf.makeH(ekf.state);
        
        %Switch back to right
        mdl.position = tmp_state_l_base(1:dof);
        mdl.velocity = tmp_state_l_base(dof+1:dof*2);
        mdl.acceleration = tmp_state_l_base(dof*2+1:end);
        mdl.forwardPosition;
        mdl.base = 'post:b_right_toe';
        %Update to latest right base results
        mdl.position = state_r(1:dof);
        mdl.velocity = state_r(dof+1:dof*2);
        mdl.acceleration = state_r(dof*2+1:end);
        mdl.forwardPosition;mdl.forwardVelocity;mdl.forwardAcceleration;
        
    end

    %Get left measures
    S_left = H_left*cov_l*H_left'+ekf.observation_noise;
    S_left = S_left(1:end-num_yaw,1:end-num_yaw);
    S_left_i = inv(S_left);
    z_myaw = z(1:end-num_yaw);
    z_hat_left_myaw = z_hat_left(1:end-num_yaw);
    eS_l(i) = (z_myaw-z_hat_left_myaw)'*S_left_i*(z_myaw-z_hat_left_myaw);
    eSN_l(i) = eS_l(i)/norm(S_left_i);
    eN_l(i) = (z_myaw-z_hat_left_myaw)'*(z_myaw-z_hat_left_myaw);
    
    %Get right measures
    S_right = H_right*cov_r*H_right'+ekf.observation_noise;
    S_right = S_right(1:end-num_yaw,1:end-num_yaw);
    S_right_i = inv(S_right);
    z_hat_right_myaw = z_hat_right(1:end-num_yaw);
    eS_r(i) = (z_myaw-z_hat_right_myaw)'*S_right_i*(z_myaw-z_hat_right_myaw);
    eSN_r(i) = eS_r(i)/norm(S_right_i);
    eN_r(i) = (z_myaw-z_hat_right_myaw)'*(z_myaw-z_hat_right_myaw);
    
    if cur_active == 'r' && eSN_r(i) - eSN_l(i) > 0 && norm(vel_l(i,:)) < 2 &&  cur_count > 20
        cur_active = 'l';
        cur_count = 0;
        %First we need to correctly set the left foot position so we use
        %the left projectes state and then we set the updated left state
        mdl.position = state_l_proj(1:dof);
        mdl.velocity = state_l_proj(dof+1:dof*2);
        mdl.acceleration = state_l_proj(dof*2+1:end);
        mdl.forwardPosition();
        mdl.base = 'post:b_left_toe';
        %Set right ankle to zero and not moving
        %indxs = [rankle_joints dof+rankle_joints dof*2+rankle_joints];
        %state_l(indxs) = 0;
        mdl.position = state_l(1:dof);
        mdl.velocity = state_l(dof+1:dof*2);
        mdl.acceleration = state_l(dof*2+1:end);
        mdl.forwardPosition();
        
        %Detach yaw sensor form the left leg and attach to swinging right
        %leg
        mdl.removeSensor(sens_yaw);
        f_yaw_base = mdl.getFrameByName('pre:b_left_upperleg:rotateZ');
        sens_yaw = SensorCore('s_yaw');
        sens_yaw.addDecorator('yaw');
        sens_yaw.base = f_yaw_base;
        f_to_b = SE3.fastInverse(f_yaw_r.t)*f_yaw_base.t;
        mdl.addSensor(sens_yaw,f_yaw_r.name,T_yaw_r);
        mdl.forwardPosition();
        
    elseif cur_active == 'l' && eSN_l(i) - eSN_r(i) > 0 && norm(vel_r(i,:)) < 2 && cur_count > 20
        cur_active = 'r';
        cur_count = 0;
        %First we need to correctly set the right foot position so we use
        %the right projectes state and then we set the updated right state
        mdl.position = state_r_proj(1:dof);
        mdl.velocity = state_r_proj(dof+1:dof*2);
        mdl.acceleration = state_r_proj(dof*2+1:end);
        mdl.forwardPosition();
        mdl.base = 'post:b_right_toe';
        %Set left ankle to zero and not moving
        %indxs = [lankle_joints dof+lankle_joints dof*2+lankle_joints];
        %state_r(indxs) = 0;
        mdl.position = state_r(1:dof);
        mdl.velocity = state_r(dof+1:dof*2);
        mdl.acceleration = state_r(dof*2+1:end);
        mdl.forwardPosition();
        
        mdl.removeSensor(sens_yaw);
        f_yaw_base = mdl.getFrameByName('pre:b_right_upperleg:rotateZ');
        sens_yaw = SensorCore('s_yaw');
        sens_yaw.addDecorator('yaw');
        sens_yaw.base = f_yaw_base;
        mdl.addSensor(sens_yaw,f_yaw_l.name,T_yaw_l);
        mdl.forwardPosition();
        
    end
    cur_count = cur_count+1;
    
    if mod(i,3) == 0 && visualize
        vis.update();
    end
    if mod(i,50) == 0 && i > 200 && visualize
        vis.addMarker(num2str(i),[foot_pos(i,:) 0]');
    end
    
    
     disp(['Frame: ' num2str(i) '/' num2str(size(mes,1)) ', Foot: ' cur_active ', el-er: ' num2str(eSN_l(i)-eSN_r(i))]);
end

figure(1);
plot(ekf_state(:, 2));


%% Look at specific peices

indx_start =5200;

for i=indx_start:size(mes,1)
    mdl.position = ekf_state(i,1:numel(mdl.joints));
    mdl.position(1:3) = 0;
    mdl.forwardPosition();
    vis.update();
    pause(dt);
end



