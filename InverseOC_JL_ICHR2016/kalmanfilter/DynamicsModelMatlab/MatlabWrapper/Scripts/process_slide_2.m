%Process the sliding along bar with one hand while other holds bar

addpath('../');
addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));

vis = rlVisualizer('vis',640,480);
vis.update();

%Load up the TRC data
bar_trc = parseTRC('C:\aslab\data\bench2\lyingdown_sliding_12_joined-bAR.trc');
bod_trc = parseTRC('C:\aslab\data\bench2\lyingdown_sliding_12_joined_smoothed-SSA_Pilot.trc');
bar_m_names = fieldnames(bar_trc.data);
bar_m_names = bar_m_names(3:end);
bod_m_names = fieldnames(bod_trc.data);
bod_m_names = bod_m_names(3:end);
%Convert to meters
for j=1:numel(bod_m_names)
    bod_trc.data.(bod_m_names{j}) = bod_trc.data.(bod_m_names{j})/1000;
end
for j=1:numel(bar_m_names)
    bar_trc.data.(bar_m_names{j})= bar_trc.data.(bar_m_names{j})/1000;
end

%Combine all trc
trc = bod_trc;
for j=1:numel(bar_m_names)
    trc.data.(bar_m_names{j})= bar_trc.data.(bar_m_names{j});
end

m_names = fieldnames( trc.data);

%% Visualize only markers
for i=1:1
    for j=1:numel(m_names)
        vis.addMarker(m_names{j},trc.data.(m_names{j})(i,:));
    end
    vis.update
end

%% Build 2 models 
mdl_rarm = rlCModel('../Models/right_arm.xml');
mdl_rarm.forwardPosition();

%Figure out the middle of the body
torso_or = (mean(bod_trc.data.SHOULDER_R)+mean(bod_trc.data.SHOULDER_L))/2;
%torso axis will point down along the body (torso origing -> clavical marker)
torso_ax = mean(bod_trc.data.CLAVICAL) - torso_or;
%Remove Z component
torso_ax = torso_ax - dot(torso_ax,[0 0 1]')*[0 0 1];
torso_ax = torso_ax/norm(torso_ax);

%Position the origin 
%NOTE WE SHIFT ALONG TORSO AXIS 4CM to get closer to SHOULDER CENTER OF
%ROTATION. Ideally we want something like harrington for shoulder.
r_shoul_or = mean(bod_trc.data.SHOULDER_R) + 0.04*torso_ax;
mdl_rarm.transforms(1).t(1:3,4) = r_shoul_or;

%Position the elbow, the primary rotation axis is determined by the vector
%between the two elbow markers lateral -> medial

%We will average this many seconds for building the model
mean_time = 1;
mean_indxs = round(1:mean_time*bod_trc.DataRate);

r_elbow_lat = mean(bod_trc.data.ELBOW_R_LAT(mean_indxs,:));
r_elbow_med = mean(bod_trc.data.ELBOW_R_MED(mean_indxs,:));
r_elbow_ax = r_elbow_med - r_elbow_lat;
r_elbow_or = r_elbow_lat + r_elbow_ax/2;
r_elbow_ax = r_elbow_ax/norm(r_elbow_ax);

r_wrist_lat = mean(bod_trc.data.WRIST_R_LAT(mean_indxs,:));
r_wrist_med = mean(bod_trc.data.WRIST_R_MED(mean_indxs,:));
r_wrist_ax = r_wrist_med - r_wrist_lat;
r_wrist_or = r_wrist_lat + r_wrist_ax/2;
r_wrist_ax = r_wrist_ax/norm(r_wrist_ax);

t_indx = find(contains({mdl_rarm.transforms.name},'length_rshoulder_relbow'),1);
trans = mdl_rarm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from shoulder to elbow in world frame
t_sho_elb_0 = r_elbow_or - r_shoul_or;
%In shoulder frame
t_sho_elb_f = f_in.t(1:3,1:3)'*t_sho_elb_0';
%Now figure out rotation between shoulder and elbow
%First figure out elbow rotation in world frame
x_ax = r_elbow_ax';
z_ax = r_wrist_or - r_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_sho_elb_f;
mdl_rarm.forwardPosition();

%Buld the elbow to wrist transform
t_indx = find(contains({mdl_rarm.transforms.name},'length_relbow_rwrist'),1);
trans = mdl_rarm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from elbow to wrist in world frame
t_elb_wri_0 = r_wrist_or - r_elbow_or;
%In shoulder frame
t_elb_wri_f = f_in.t(1:3,1:3)'*t_elb_wri_0';

%Now figure out rotation between elbow and wrist, note, we keep same Z axis
%as before going from the elbow up to wrist and through
x_ax = r_wrist_ax';
z_ax = r_wrist_or - r_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_elb_wri_f;
mdl_rarm.forwardPosition();

%Finally make the end effector the point on the bar closest to the IMU

%Figure out the bar axis
right_u = mean(bar_trc.data.RIGHT_U(mean_indxs,:));
right_l = mean(bar_trc.data.RIGHT_L(mean_indxs,:));
r_bar_end = (right_u+right_l)/2;
left_u = mean(bar_trc.data.LEFT_U(mean_indxs,:));
left_l = mean(bar_trc.data.LEFT_L(mean_indxs,:));
l_bar_end = (left_u+left_l)/2;
bar_vec = r_bar_end-l_bar_end;
bar_axis = bar_vec/norm(bar_vec);
%So any point on the bar is expressed as l_bar_end + t*bar_axis
%Lets project the right hand IMU onto the bar 
r_wrist_imu = mean(trc.data.WRIST_R_IMU_T(mean_indxs,:) + ...
    trc.data.WRIST_R_IMU_L(mean_indxs,:) + ...
    trc.data.WRIST_R_IMU_R(mean_indxs,:))/3;
l_bar_end_to_r_hand = r_wrist_imu - l_bar_end;
a = dot(l_bar_end_to_r_hand,bar_axis);
r_hand = l_bar_end + a*bar_axis;

%Make the right hand Y axis aligned with the bar
%Rotation from left hand to world
R0_rh = eye(3);
R0_rh(:,2) = bar_axis;
%Z = x cross y
R0_rh(:,3) = cross([1 0 0]',R0_rh(:,2));
R0_rh(:,3) = R0_rh(:,3)/norm(R0_rh(:,3));
%X = y cross z 
R0_rh(:,1) = cross(R0_rh(:,2),R0_rh(:,3));
R0_rh(:,1) = R0_rh(:,1)/norm(R0_rh(:,1));

t_indx = find(contains({mdl_rarm.transforms.name},'length_rwrist_rhand'),1);
trans = mdl_rarm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

t_wri_r_bar_end = r_hand - r_wrist_or;
t_wri_bar_f = f_in.t(1:3,1:3)'*t_wri_r_bar_end';
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R0_lh;
trans.t(1:3,4) = t_wri_bar_f;
mdl_rarm.forwardPosition();

%Now in the same exact way build the left arm
mdl_larm = rlCModel('../Models/left_arm.xml');
mdl_larm.forwardPosition();

%Position the origin
%SIMILARLY SHIFT A BIT ALONG TORSO AXIS
l_shoul_or = mean(bod_trc.data.SHOULDER_L) + 0.04*torso_ax;
mdl_larm.transforms(1).t(1:3,4) = l_shoul_or;

%Position the elbow, the primary rotation axis is determined by the vector
%between the two elbow markers lateral -> medial

%We will average this many seconds for building the model
mean_time = 2;
mean_indxs = round(1:mean_time*bod_trc.DataRate);

l_elbow_lat = mean(bod_trc.data.ELBOW_L_LAT(mean_indxs,:));
l_elbow_med = mean(bod_trc.data.ELBOW_L_MED(mean_indxs,:));
l_elbow_ax = l_elbow_med - l_elbow_lat;
l_elbow_or = l_elbow_lat + l_elbow_ax/2;
l_elbow_ax = l_elbow_ax/norm(l_elbow_ax);

l_wrist_lat = mean(bod_trc.data.WRIST_L_LAT(mean_indxs,:));
l_wrist_med = mean(bod_trc.data.WRIST_L_MED(mean_indxs,:));
l_wrist_ax = l_wrist_med - l_wrist_lat;
l_wrist_or = l_wrist_lat + l_wrist_ax/2;
l_wrist_ax = l_wrist_ax/norm(l_wrist_ax);

t_indx = find(contains({mdl_larm.transforms.name},'length_lshoulder_lelbow'),1);
trans = mdl_larm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from shoulder to elbow in world frame
t_sho_elb_0 = l_elbow_or - l_shoul_or;
%In shoulder frame
t_sho_elb_f = f_in.t(1:3,1:3)'*t_sho_elb_0';
%Now figure out rotation between shoulder and elbow
%First figure out elbow rotation in world frame
x_ax = l_elbow_ax';
z_ax = l_wrist_or - l_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_sho_elb_f;
mdl_larm.forwardPosition();

%Buld the elbow to wrist transform
t_indx = find(contains({mdl_larm.transforms.name},'length_lelbow_lwrist'),1);
trans = mdl_larm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from elbow to wrist in world frame
t_elb_wri_0 = l_wrist_or - l_elbow_or;
%In shoulder frame
t_elb_wri_f = f_in.t(1:3,1:3)'*t_elb_wri_0';

%Now figure out rotation between elbow and wrist, note, we keep same Z axis
%as before going from the elbow up to wrist and through
x_ax = l_wrist_ax';
z_ax = l_wrist_or - l_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_elb_wri_f;
mdl_larm.forwardPosition();

l_wrist_imu = mean(trc.data.WRIST_L_IMU_T(mean_indxs,:) + ...
    trc.data.WRIST_L_IMU_L(mean_indxs,:) + ...
    trc.data.WRIST_L_IMU_R(mean_indxs,:))/3;
l_bar_end_to_l_hand = l_wrist_imu - l_bar_end;
a = dot(l_bar_end_to_l_hand,bar_axis);
l_hand = l_bar_end + a*bar_axis;

%Make the left hand Y axis aligned with the bar
%Rotation from left hand to world
R0_lh = eye(3);
R0_lh(:,2) = bar_axis;
%Z = x cross y
R0_lh(:,3) = cross([1 0 0]',R0_lh(:,2));
R0_lh(:,3) = R0_lh(:,3)/norm(R0_lh(:,3));
%X = y cross z 
R0_lh(:,1) = cross(R0_lh(:,2),R0_lh(:,3));
R0_lh(:,1) = R0_lh(:,1)/norm(R0_lh(:,1));

t_indx = find(contains({mdl_larm.transforms.name},'length_lwrist_lhand'),1);
trans = mdl_larm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

t_wri_l_bar_end = l_hand - l_wrist_or;
t_wri_bar_f = f_in.t(1:3,1:3)'*t_wri_l_bar_end';
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R0_lh;
trans.t(1:3,4) = t_wri_bar_f;
mdl_larm.forwardPosition();

%% Attach Motion Capture Markers

%This cell array defines model frames and what markers are attached to them
marker_attachement_rarm = ...
    {{'body_rshoulder_relbow','ELBOW_R_LAT','ELBOW_R_MED','UPPERARM_R_IMU_T','UPPERARM_R_IMU_L','UPPERARM_R_IMU_R'};...
    {'body_relbow_rwrist','WRIST_R_LAT','WRIST_R_MED','LOWERARM_R_IMU_T','LOWERARM_R_IMU_R','LOWERARM_R_IMU_L'};...
    ...%{'body_rwrist_rhand','RIGHT_U','RIGHT_L'};
    };

for i=1:numel(marker_attachement_rarm)
    f_name =  marker_attachement_rarm{i}{1};
    frame = mdl_rarm.getFrameByName(f_name);
    
    for j = 2:numel(marker_attachement_rarm{i})
        %translation from frame to marker in world frame
        marker = mean(trc.data.(marker_attachement_rarm{i}{j})(mean_indxs,:))';
        t_0 = marker - frame.t(1:3,4);
        t_f = frame.t(1:3,1:3)'*t_0;
        T = eye(4);
        T(1:3,4) = t_f;
        sens = SensorCore(marker_attachement_rarm{i}{j});
        sens.addDecorator('position');
        mdl_rarm.addSensor(sens,frame.name,T);
    end
end

mdl_rarm.forwardPosition();

marker_attachement_larm = ...
    {{'body_lshoulder_lelbow','ELBOW_L_LAT','ELBOW_L_MED','UPPERARM_L_IMU_T','UPPERARM_L_IMU_L','UPPERARM_L_IMU_R'};...
    {'body_lelbow_lwrist','WRIST_L_LAT','WRIST_L_MED','LOWERARM_L_IMU_T','LOWERARM_L_IMU_R','LOWERARM_L_IMU_L'};...
    {'body_lwrist_lhand','LEFT_U','LEFT_L'};
    };

for i=1:numel(marker_attachement_larm)
    f_name =  marker_attachement_larm{i}{1};
    frame = mdl_larm.getFrameByName(f_name);
    
    for j = 2:numel(marker_attachement_larm{i})
        %translation from frame to marker in world frame
        marker = mean(trc.data.(marker_attachement_larm{i}{j})(mean_indxs,:))';
        t_0 = marker - frame.t(1:3,4);
        t_f = frame.t(1:3,1:3)'*t_0;
        T = eye(4);
        T(1:3,4) = t_f;
        sens = SensorCore(marker_attachement_larm{i}{j});
        sens.addDecorator('position');
        mdl_larm.addSensor(sens,frame.name,T);
    end
end
mdl_larm.forwardPosition();

%Build measurement vectors
mes_rarm = [];
for i=1:numel(marker_attachement_rarm)
    for j = 2:numel(marker_attachement_rarm{i})
        mes_rarm = [mes_rarm trc.data.(marker_attachement_rarm{i}{j})];
    end
end
mes_larm = [];
for i=1:numel(marker_attachement_larm)
    for j = 2:numel(marker_attachement_larm{i})
        mes_larm = [mes_larm trc.data.(marker_attachement_larm{i}{j})];
    end
end

mes_rarm_obj = SensorMeasurement(trc.NumFrames,numel(mdl_rarm.sensors));
[mes_rarm_obj.type] = deal(1);
[mes_rarm_obj.size] = deal(3);
mes_rarm_obj.setMesArray(mes_rarm);
% marker_attachement_larm
mes_larm_obj = SensorMeasurement(trc.NumFrames,numel(mdl_larm.sensors));
[mes_larm_obj.type] = deal(1);
[mes_larm_obj.size] = deal(3);
mes_larm_obj.setMesArray(mes_larm);

vis.addModel(mdl_rarm);
vis.addModel(mdl_larm);
vis.update();

%% Set up EKF to estimate joint angles based on MOCAP
ekf_rarm = EKF_Q_DQ_DDQ(mdl_rarm);
ekf_larm = EKF_Q_DQ_DDQ(mdl_larm);
state_est_rarm_mocap = zeros(trc.NumFrames,numel(ekf_rarm.state));
state_est_larm_mocap = zeros(trc.NumFrames,numel(ekf_larm.state));
%Marker Noise
ekf_rarm.observation_noise = ekf_rarm.observation_noise*0.01;
ekf_larm.observation_noise = ekf_larm.observation_noise*0.01;

eta = 10;
dt = 1/trc.DataRate;
%Process noise
dim = numel(mdl_rarm.joints);
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
ekf_larm.process_noise = P;

%% Do the estimation

for i=1:trc.NumFrames
    
    %This is the measurement for current timestep, must be vertical vector
    z_rarm = mes_rarm_obj(i,:);
    z_larm = mes_larm_obj(i,:);
    
    ekf_rarm.run_iteration(dt,z_rarm);
    state_est_rarm_mocap(i,:) = ekf_rarm.state;
    
    ekf_larm.run_iteration(dt,z_larm);
    state_est_larm_mocap(i,:) = ekf_larm.state;
    
    %Draw markers
    for j=1:numel(m_names)
        vis.addMarker(m_names{j},trc.data.(m_names{j})(i,:));
    end
    vis.update
end

%% Now same using IMU data 

%% Build 2 models 
mdl_rarm = rlCModel('../Models/right_arm.xml');
mdl_rarm.forwardPosition();

%Figure out the middle of the body
torso_or = (mean(bod_trc.data.SHOULDER_R)+mean(bod_trc.data.SHOULDER_L))/2;
%torso axis will point down along the body (torso origing -> clavical marker)
torso_ax = mean(bod_trc.data.CLAVICAL) - torso_or;
%Remove Z component
torso_ax = torso_ax - dot(torso_ax,[0 0 1]')*[0 0 1];
torso_ax = torso_ax/norm(torso_ax);

%Position the origin 
%NOTE WE SHIFT ALONG TORSO AXIS 4CM to get closer to SHOULDER CENTER OF
%ROTATION. Ideally we want something like harrington for shoulder.
r_shoul_or = mean(bod_trc.data.SHOULDER_R) + 0.04*torso_ax;
mdl_rarm.transforms(1).t(1:3,4) = r_shoul_or;

%Position the elbow, the primary rotation axis is determined by the vector
%between the two elbow markers lateral -> medial

%We will average this many seconds for building the model
mean_time = 1;
mean_indxs = round(1:mean_time*bod_trc.DataRate);

r_elbow_lat = mean(bod_trc.data.ELBOW_R_LAT(mean_indxs,:));
r_elbow_med = mean(bod_trc.data.ELBOW_R_MED(mean_indxs,:));
r_elbow_ax = r_elbow_med - r_elbow_lat;
r_elbow_or = r_elbow_lat + r_elbow_ax/2;
r_elbow_ax = r_elbow_ax/norm(r_elbow_ax);

r_wrist_lat = mean(bod_trc.data.WRIST_R_LAT(mean_indxs,:));
r_wrist_med = mean(bod_trc.data.WRIST_R_MED(mean_indxs,:));
r_wrist_ax = r_wrist_med - r_wrist_lat;
r_wrist_or = r_wrist_lat + r_wrist_ax/2;
r_wrist_ax = r_wrist_ax/norm(r_wrist_ax);

t_indx = find(contains({mdl_rarm.transforms.name},'length_rshoulder_relbow'),1);
trans = mdl_rarm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from shoulder to elbow in world frame
t_sho_elb_0 = r_elbow_or - r_shoul_or;
%In shoulder frame
t_sho_elb_f = f_in.t(1:3,1:3)'*t_sho_elb_0';
%Now figure out rotation between shoulder and elbow
%First figure out elbow rotation in world frame
x_ax = r_elbow_ax';
z_ax = r_wrist_or - r_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_sho_elb_f;
mdl_rarm.forwardPosition();

%Buld the elbow to wrist transform
t_indx = find(contains({mdl_rarm.transforms.name},'length_relbow_rwrist'),1);
trans = mdl_rarm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from elbow to wrist in world frame
t_elb_wri_0 = r_wrist_or - r_elbow_or;
%In shoulder frame
t_elb_wri_f = f_in.t(1:3,1:3)'*t_elb_wri_0';

%Now figure out rotation between elbow and wrist, note, we keep same Z axis
%as before going from the elbow up to wrist and through
x_ax = r_wrist_ax';
z_ax = r_wrist_or - r_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_elb_wri_f;
mdl_rarm.forwardPosition();

%Finally make the end effector the point on the bar closest to the IMU

%Figure out the bar axis
right_u = mean(bar_trc.data.RIGHT_U(mean_indxs,:));
right_l = mean(bar_trc.data.RIGHT_L(mean_indxs,:));
r_bar_end = (right_u+right_l)/2;
left_u = mean(bar_trc.data.LEFT_U(mean_indxs,:));
left_l = mean(bar_trc.data.LEFT_L(mean_indxs,:));
l_bar_end = (left_u+left_l)/2;
bar_vec = r_bar_end-l_bar_end;
bar_axis = bar_vec/norm(bar_vec);
%So any point on the bar is expressed as l_bar_end + t*bar_axis
%Lets project the right hand IMU onto the bar 
r_wrist_imu = mean(trc.data.WRIST_R_IMU_T(mean_indxs,:) + ...
    trc.data.WRIST_R_IMU_L(mean_indxs,:) + ...
    trc.data.WRIST_R_IMU_R(mean_indxs,:))/3;
l_bar_end_to_r_hand = r_wrist_imu - l_bar_end;
a = dot(l_bar_end_to_r_hand,bar_axis);
r_hand = l_bar_end + a*bar_axis;

t_indx = find(contains({mdl_rarm.transforms.name},'length_rwrist_rhand'),1);
trans = mdl_rarm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

t_wri_r_bar_end = r_hand - r_wrist_or;
t_wri_bar_f = f_in.t(1:3,1:3)'*t_wri_r_bar_end';
trans.t(1:3,1:3) = f_in.t(1:3,1:3)';
trans.t(1:3,4) = t_wri_bar_f;
mdl_rarm.forwardPosition();

%Now in the same exact way build the left arm
mdl_larm = rlCModel('../Models/left_arm_bar.xml');
mdl_larm.forwardPosition();

%Position the origin
%SIMILARLY SHIFT A BIT ALONG TORSO AXIS
l_shoul_or = mean(bod_trc.data.SHOULDER_L) + 0.04*torso_ax;
mdl_larm.transforms(1).t(1:3,4) = l_shoul_or;

%Position the elbow, the primary rotation axis is determined by the vector
%between the two elbow markers lateral -> medial

%We will average this many seconds for building the model
mean_time = 2;
mean_indxs = round(1:mean_time*bod_trc.DataRate);

l_elbow_lat = mean(bod_trc.data.ELBOW_L_LAT(mean_indxs,:));
l_elbow_med = mean(bod_trc.data.ELBOW_L_MED(mean_indxs,:));
l_elbow_ax = l_elbow_med - l_elbow_lat;
l_elbow_or = l_elbow_lat + l_elbow_ax/2;
l_elbow_ax = l_elbow_ax/norm(l_elbow_ax);

l_wrist_lat = mean(bod_trc.data.WRIST_L_LAT(mean_indxs,:));
l_wrist_med = mean(bod_trc.data.WRIST_L_MED(mean_indxs,:));
l_wrist_ax = l_wrist_med - l_wrist_lat;
l_wrist_or = l_wrist_lat + l_wrist_ax/2;
l_wrist_ax = l_wrist_ax/norm(l_wrist_ax);

t_indx = find(contains({mdl_larm.transforms.name},'length_lshoulder_lelbow'),1);
trans = mdl_larm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from shoulder to elbow in world frame
t_sho_elb_0 = l_elbow_or - l_shoul_or;
%In shoulder frame
t_sho_elb_f = f_in.t(1:3,1:3)'*t_sho_elb_0';
%Now figure out rotation between shoulder and elbow
%First figure out elbow rotation in world frame
x_ax = l_elbow_ax';
z_ax = l_wrist_or - l_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_sho_elb_f;
mdl_larm.forwardPosition();

%Buld the elbow to wrist transform
t_indx = find(contains({mdl_larm.transforms.name},'length_lelbow_lwrist'),1);
trans = mdl_larm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

%Translation from elbow to wrist in world frame
t_elb_wri_0 = l_wrist_or - l_elbow_or;
%In shoulder frame
t_elb_wri_f = f_in.t(1:3,1:3)'*t_elb_wri_0';

%Now figure out rotation between elbow and wrist, note, we keep same Z axis
%as before going from the elbow up to wrist and through
x_ax = l_wrist_ax';
z_ax = l_wrist_or - l_elbow_or;
z_ax = z_ax'/norm(z_ax);
y_ax = cross(z_ax,x_ax);
y_ax = y_ax/norm(y_ax);
z_ax = cross(x_ax,y_ax);
R = [x_ax y_ax z_ax];
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R;
trans.t(1:3,4) = t_elb_wri_f;
mdl_larm.forwardPosition();

l_wrist_imu = mean(trc.data.WRIST_L_IMU_T(mean_indxs,:) + ...
    trc.data.WRIST_L_IMU_L(mean_indxs,:) + ...
    trc.data.WRIST_L_IMU_R(mean_indxs,:))/3;
l_bar_end_to_l_hand = l_wrist_imu - l_bar_end;
a = dot(l_bar_end_to_l_hand,bar_axis);
l_hand = l_bar_end + a*bar_axis;

%Make the left hand Y axis aligned with the bar
%Rotation from left hand to world
R0_lh = eye(3);
R0_lh(:,2) = bar_axis;
%Z = x cross y
R0_lh(:,3) = cross([1 0 0]',R0_lh(:,2));
R0_lh(:,3) = R0_lh(:,3)/norm(R0_lh(:,3));
%X = y cross z 
R0_lh(:,1) = cross(R0_lh(:,2),R0_lh(:,3));
R0_lh(:,1) = R0_lh(:,1)/norm(R0_lh(:,1));

t_indx = find(contains({mdl_larm.transforms.name},'length_lwrist_lhand'),1);
trans = mdl_larm.transforms(t_indx);
f_in = trans.frame_in;
f_out = trans.frame_out;

t_wri_l_bar_end = l_hand - l_wrist_or;
t_wri_bar_f = f_in.t(1:3,1:3)'*t_wri_l_bar_end';
trans.t(1:3,1:3) = f_in.t(1:3,1:3)'*R0_lh;
trans.t(1:3,4) = t_wri_bar_f;
mdl_larm.forwardPosition();

%Add the ends of the bar on each side 
t_indx = find(contains({mdl_larm.transforms.name},'length_lhand_lbar'),1);
trans = mdl_larm.transforms(t_indx);
%We know that the bar is aligned along Y frame of the lhand so we just need
%the distance 
l_hand_to_lbarend = -norm(l_bar_end-l_hand);
trans.t(2,4) = l_hand_to_lbarend;
%Same for right
t_indx = find(contains({mdl_larm.transforms.name},'length_lhand_rbar'),1);
trans = mdl_larm.transforms(t_indx);
%We know that the bar is aligned along Y frame of the lhand so we just need
%the distance 
l_hand_to_rbarend = norm(r_bar_end-l_hand);
trans.t(2,4) = l_hand_to_rbarend;

%% Attach the IMUS 

%This defines the IMU attachement {frame to attach to, imu markers}
imu_attachement_larm = {...
    {'body_lshoulder_lelbow','UPPERARM_L_IMU_L','UPPERARM_L_IMU_T','UPPERARM_L_IMU_R'};...
    {'body_lelbow_lwrist','LOWERARM_L_IMU_L','LOWERARM_L_IMU_T','LOWERARM_L_IMU_R'};...
    {'body_lwrist_lhand','WRIST_L_IMU_L','WRIST_L_IMU_T','WRIST_L_IMU_R'};};

for i=1:numel(imu_attachement_larm)
    f_name =  imu_attachement_larm{i}{1};
    frame = mdl_larm.getFrameByName(f_name);
    
    P = mean(trc.data.(imu_attachement_larm{i}{2})(mean_indxs,:));
    Q = mean(trc.data.(imu_attachement_larm{i}{3})(mean_indxs,:));
    R = mean(trc.data.(imu_attachement_larm{i}{4})(mean_indxs,:));
    [R_imu_0, R_0_imu] = points2rot(P,Q,R);
    t_0 =  mean([P;Q;R])' - frame.t(1:3,4);
    t_f = frame.t(1:3,1:3)'*t_0;
    
    sens = SensorCore(imu_attachement_larm{i}{2}(1:end-2));
    sens.addDecorator('gyroscope');
    sens.addDecorator('accelerometer');
    T = eye(4);
    T(1:3,1:3) = frame.t(1:3,1:3)'*R_0_imu';
    T(1:3,4) = t_f;
    mdl_larm.addSensor(sens,frame.name,T)
end
mdl_larm.forwardPosition();


imu_attachement_rarm = {...
    {'body_rshoulder_relbow','UPPERARM_R_IMU_L','UPPERARM_R_IMU_T','UPPERARM_R_IMU_R'};...
    {'body_relbow_rwrist','LOWERARM_R_IMU_L','LOWERARM_R_IMU_T','LOWERARM_R_IMU_R'};...
    {'body_rwrist_rhand','WRIST_R_IMU_L','WRIST_R_IMU_T','WRIST_R_IMU_R'};...
    };

for i=1:numel(imu_attachement_rarm)
    f_name =  imu_attachement_rarm{i}{1};
    frame = mdl_rarm.getFrameByName(f_name);
    
    P = mean(trc.data.(imu_attachement_rarm{i}{2})(mean_indxs,:));
    Q = mean(trc.data.(imu_attachement_rarm{i}{3})(mean_indxs,:));
    R = mean(trc.data.(imu_attachement_rarm{i}{4})(mean_indxs,:));
    [R_imu_0, R_0_imu] = points2rot(P,Q,R);
    t_0 =  mean([P;Q;R])' - frame.t(1:3,4);
    t_f = frame.t(1:3,1:3)'*t_0;
    
    sens = SensorCore(imu_attachement_rarm{i}{2}(1:end-2));
    sens.addDecorator('gyroscope');
    sens.addDecorator('accelerometer');
    T = eye(4);
    T(1:3,1:3) = frame.t(1:3,1:3)'*R_0_imu';
    T(1:3,4) = t_f;
    mdl_rarm.addSensor(sens,frame.name,T)
end
mdl_rarm.forwardPosition();

%% Set upe IMU measurement Vectors 
time_stamp = '2019_02_26_12_17_48';
imu_dat_path = 'C:\aslab\data\bench2\imu\';
r_arm_imu_files = {['0006667D7185_' time_stamp '.header'],...
    ['00066681EB52_' time_stamp '.header'],...
    ['0006667D713B_' time_stamp '.header'],...
    };
imus_rarm = tinySensorDataHandle.empty();
num_samples_rarm = [];
for i=1:numel(r_arm_imu_files)
    imu = arsLoader([imu_dat_path r_arm_imu_files{i}]);
    imus_rarm(i) = imu;
    num_samples_rarm(i) = size(imu.accelerometerCalibrated,1);
end

l_arm_imu_files = {['00066681EA9C_' time_stamp '.header'],...
    ['00066681EA9B_' time_stamp '.header'],...
    ['00066681EA98_' time_stamp '.header']};
imus_larm = tinySensorDataHandle.empty();
num_samples_larm = [];
for i=1:numel(l_arm_imu_files)
    imu = arsLoader([imu_dat_path l_arm_imu_files{i}]);
    imus_larm(i) = imu;
    num_samples_larm(i) = size(imu.accelerometerCalibrated,1);
end

min_samples = min([num_samples_rarm num_samples_larm]);

mes_larm = [];
for i = 1:numel(imus_larm)
    mes_larm = [mes_larm imus_larm(i).gyroscopeCalibrated(1:min_samples,:)...
        imus_larm(i).accelerometerCalibrated(1:min_samples,:)];
end
mes_rarm = [];
for i = 1:numel(imus_rarm)
    mes_rarm = [mes_rarm imus_rarm(i).gyroscopeCalibrated(1:min_samples,:)...
        imus_rarm(i).accelerometerCalibrated(1:min_samples,:)];
end

%% Set up and run Multi Chain EKF 

%Reset model
mdl_rarm.position(:) = 0;mdl_rarm.velocity(:) = 0;mdl_rarm.acceleration(:) = 0;
mdl_larm.position(:) = 0;mdl_larm.velocity(:) = 0;mdl_larm.acceleration(:) = 0;
mdl_rarm.forwardPosition();
mdl_larm.forwardPosition();

ekf = MC_EKF_Q_DQ_DDQ();
ekf.addModel(mdl_larm);
ekf.addModel(mdl_rarm);
ekf.observation_noise = diag(repmat([0.1 0.1 0.1 10 10 10],1,6));
ekf.covariance = ekf.process_noise;
state_est_imu = zeros(size(mes_rarm,1),ekf.sizeX);

%Add Constraints between the two chains
C1 = EKF_Constraint_T0EE('C1',mdl_larm,'frame_lhand_end',eye(4),mdl_rarm,'frame_rhand_end',eye(4));
C1.type(:) = false;C1.type = logical([1 0 0 1 0 1]');
ekf.addConstraint(C1);

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl_rarm);
vis.addModel(mdl_larm);
%Set Constraint Colors of constraints
const_color = [0 1 1 1]';
for i=1:numel(ekf.constraints)
    C = ekf.constraints(i);
    vis.setSensorColor(C.models(1),C.sensors(1),const_color);
    vis.setSensorColor(C.models(2),C.sensors(2),const_color);
end
%Set the bar color
t_indx = find(contains({mdl_larm.transforms.name},'length_lhand_rbar'),1);
trans = mdl_larm.transforms(t_indx);
vis.setTransformColor(mdl_larm,trans,[0.5 0.5 0.5 1]');
t_indx = find(contains({mdl_larm.transforms.name},'length_lhand_lbar'),1);
trans = mdl_larm.transforms(t_indx);
vis.setTransformColor(mdl_larm,trans,[0.5 0.5 0.5 1]');
%Set visualizer background
vis.setBackgroundColour([1 1 1]);

vis.update();

%% Test the constraint Jacobian computatino 
mdl_larm.forwardPosition();
mdl_rarm.forwardPosition();
f_l = mdl_larm.getFrameByName('frame_lhand_end');
f_r = mdl_rarm.getFrameByName('frame_rhand_end');

J_num = zeros(6,numel(mdl_larm.joints) + numel(mdl_rarm.joints));
dT0_eel_num = zeros(4,4,numel(mdl_larm.joints));
dTeel_0_num = zeros(4,4,numel(mdl_larm.joints));
dT0_eer_num = zeros(4,4,numel(mdl_larm.joints));
dTeer_0_num = zeros(4,4,numel(mdl_larm.joints));

ep = 1e-6;
for i=1:numel(mdl_larm.joints)
    
    %Transform from right end effector to left
    T0_eel = f_l.t;
    Teel_0 = SE3.fastInverse(T0_eel);
    T0_eer = f_r.t;
    Teer_0 = SE3.fastInverse(T0_eer);
    
    Teel_eer = SE3.fastInverse(f_l.t)*f_r.t;
    pos_init = mdl_larm.position(i);
    mdl_larm.position(i) = mdl_larm.position(i) + ep;
    mdl_larm.forwardPosition();
    T0_eel_eps = f_l.t;
    Teel_0_eps = SE3.fastInverse(T0_eel_eps);
    T0_eer_eps = f_r.t;
    Teer_0_eps = SE3.fastInverse(T0_eer_eps);
    
    Teel_eer_eps = SE3.fastInverse(f_l.t)*f_r.t;
    
    J_num(1:3,i) = (T0_eel_eps(1:3,C1.type(1:3)) - T0_eel(1:3,C1.type(1:3)))/ep;
    J_num(4:6,i) = (Teel_eer_eps(1:3,4) - Teel_eer(1:3,4))/ep;
    dT0_eel_num(:,:,i) = (T0_eel_eps - T0_eel)/ep;
    dTeel_0_num(:,:,i) = (Teel_0_eps - Teel_0)/ep;
    
    mdl_larm.position(i)= pos_init;
    mdl_larm.forwardPosition();
    
end


%% Run EKF and see what happens
dt = 1/30;

f_l = mdl_larm.getFrameByName('frame_lhand_end');
f_r = mdl_rarm.getFrameByName('frame_rhand_end');

%Transformation from right EE to left EE in left EE frame
Tl_r = zeros(4,4,size(mes_rarm,1));

for i=1:size(mes_rarm,1)
    ekf.run_iteration(dt,[mes_larm(i,:) mes_rarm(i,:)]');
    state_est_imu(i,:) = ekf.state;
    
    Tl_r(:,:,i) = SE3.fastInverse(f_l.t)*f_r.t;
    
    vis.update();
end

%% Lets look at plots

%Translation of right end effector in left end effector frame
tmp = permute(Tl_r,[1,3,2]);
t_l_r = tmp(1:3,:,4)';

%X axis of EE 2 in frame of EE 1
x_l_r = tmp(1:3,:,1)';


