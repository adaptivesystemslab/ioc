%Process the sliding along bar with one hand while other holds bar

addpath('../');
addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));

%Load post alignment data
load('C:\aslab\data\bench2\featureSet_mocap_post_const1_Subject02_lyingdown_sliding_12.mat')
trc = dataInstance_trc_post;
m_names = fieldnames( trc.data);

%% Build 2 models 
mdl_rarm = rlCModel('../Models/right_arm.xml');
mdl_rarm.forwardPosition();

%Figure out the middle of the body
torso_or = (mean(trc.data.SHOULDER_R)+mean(trc.data.SHOULDER_L))/2;
%torso axis will point down along the body (torso origing -> clavical marker)
torso_ax = mean(trc.data.CLAVICAL) - torso_or;
%Remove Z component
torso_ax = torso_ax - dot(torso_ax,[0 0 1]')*[0 0 1];
torso_ax = torso_ax/norm(torso_ax);

%Position the origin 
%NOTE WE SHIFT ALONG TORSO AXIS 4CM to get closer to SHOULDER CENTER OF
%ROTATION. Ideally we want something like harrington for shoulder.
r_shoul_or = mean(trc.data.SHOULDER_R) + 0.04*torso_ax;
mdl_rarm.transforms(1).t(1:3,4) = r_shoul_or;

%Position the elbow, the primary rotation axis is determined by the vector
%between the two elbow markers lateral -> medial

%We will average this many seconds for building the model
mean_time = 1;
mean_indxs = round(1:mean_time*trc.DataRate);

r_elbow_lat = mean(trc.data.ELBOW_R_LAT(mean_indxs,:));
r_elbow_med = mean(trc.data.ELBOW_R_MED(mean_indxs,:));
r_elbow_ax = r_elbow_med - r_elbow_lat;
r_elbow_or = r_elbow_lat + r_elbow_ax/2;
r_elbow_ax = r_elbow_ax/norm(r_elbow_ax);

r_wrist_lat = mean(trc.data.WRIST_R_LAT(mean_indxs,:));
r_wrist_med = mean(trc.data.WRIST_R_MED(mean_indxs,:));
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
right_u = mean(trc.data.RIGHT_U(mean_indxs,:));
right_l = mean(trc.data.RIGHT_L(mean_indxs,:));
r_bar_end = (right_u+right_l)/2;
left_u = mean(trc.data.LEFT_U(mean_indxs,:));
left_l = mean(trc.data.LEFT_L(mean_indxs,:));
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
l_shoul_or = mean(trc.data.SHOULDER_L) + 0.04*torso_ax;
mdl_larm.transforms(1).t(1:3,4) = l_shoul_or;

%Position the elbow, the primary rotation axis is determined by the vector
%between the two elbow markers lateral -> medial

%We will average this many seconds for building the model
mean_time = 2;
mean_indxs = round(1:mean_time*trc.DataRate);

l_elbow_lat = mean(trc.data.ELBOW_L_LAT(mean_indxs,:));
l_elbow_med = mean(trc.data.ELBOW_L_MED(mean_indxs,:));
l_elbow_ax = l_elbow_med - l_elbow_lat;
l_elbow_or = l_elbow_lat + l_elbow_ax/2;
l_elbow_ax = l_elbow_ax/norm(l_elbow_ax);

l_wrist_lat = mean(trc.data.WRIST_L_LAT(mean_indxs,:));
l_wrist_med = mean(trc.data.WRIST_L_MED(mean_indxs,:));
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

%% Visualize the model without IMUs and fade them in
vis = rlVisualizer('vis',1920,1080);
vis.addModel(mdl_rarm);
vis.addModel(mdl_larm);
vis.setBackgroundColour([1 1 1]);

viewport_init = [-1.0860    0.1624         0   -2.8000   13.4001    3.5274]';
vis.setViewportPose(viewport_init);

vis.update();

for j=1:numel(mdl_rarm.sensors)
    vis.setSensorColor(mdl_rarm,mdl_rarm.sensors(j),[1 0 1 0]);
    t_indx = find(contains({mdl_rarm.transforms.name},[mdl_rarm.sensors(j).name 't']),1);
    trans = mdl_rarm.transforms(t_indx);
    vis.setTransformColor(mdl_rarm,trans,[0 0 1 0]);
end
for j=1:numel(mdl_larm.sensors)
    vis.setSensorColor(mdl_larm,mdl_larm.sensors(j),[1 0 1 0]);
    t_indx = find(contains({mdl_larm.transforms.name},[mdl_larm.sensors(j).name 't']),1);
    trans = mdl_larm.transforms(t_indx);
    vis.setTransformColor(mdl_larm,trans,[0 0 1 0]);
end


%% Fade in BAR
t_indx = find(contains({mdl_larm.transforms.name},'length_lhand_rbar'),1);
trans_r = mdl_larm.transforms(t_indx);
vis.setTransformColor(mdl_larm,trans_r,[0.5 0.5 0.5 0]');
t_indx = find(contains({mdl_larm.transforms.name},'length_lhand_lbar'),1);
trans_l = mdl_larm.transforms(t_indx);
vis.setTransformColor(mdl_larm,trans_l,[0.5 0.5 0.5 0]');

vis.update();

vid_path = 'C:\asl_git\kalmanfilter\General_FKEKF\Paper_Const\Video\';
vid_file = 'fade_bar.avi';

vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=0:0.005:1
    vis.setTransformColor(mdl_larm,trans_r,[0.5 0.5 0.5 i]);
    vis.setTransformColor(mdl_larm,trans_l,[0.5 0.5 0.5 i]);
    vis.update();
    b = false;
    while(~b)
        [I, b] = vis.getScreenshot();
    end
    writeVideo(vid,I);
end
close(vid);

%% Fade in IMUS

vid_file = 'fade_imu.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for j=1:numel(mdl_rarm.sensors)
    vis.setSensorColor(mdl_rarm,mdl_rarm.sensors(j),[1 0 1 0]);
    t_indx = find(contains({mdl_rarm.transforms.name},[mdl_rarm.sensors(j).name 't']),1);
    trans = mdl_rarm.transforms(t_indx);
    vis.setTransformColor(mdl_rarm,trans,[0 0 1 0]);
end
for j=1:numel(mdl_larm.sensors)
    vis.setSensorColor(mdl_larm,mdl_larm.sensors(j),[1 0 1 0]);
    t_indx = find(contains({mdl_larm.transforms.name},[mdl_larm.sensors(j).name 't']),1);
    trans = mdl_larm.transforms(t_indx);
    vis.setTransformColor(mdl_larm,trans,[0 0 1 0]);
end
vis.update();

for i=0:0.005:1
    for j=1:numel(mdl_rarm.sensors)
       vis.setSensorColor(mdl_rarm,mdl_rarm.sensors(j),[1 0 1 i]); 
       t_indx = find(contains({mdl_rarm.transforms.name},[mdl_rarm.sensors(j).name 't']),1);
       trans = mdl_rarm.transforms(t_indx);
       vis.setTransformColor(mdl_rarm,trans,[0 0 1 i]); 
    end
    for j=1:numel(mdl_larm.sensors)
       vis.setSensorColor(mdl_larm,mdl_larm.sensors(j),[1 0 1 i]); 
       t_indx = find(contains({mdl_larm.transforms.name},[mdl_larm.sensors(j).name 't']),1);
       trans = mdl_larm.transforms(t_indx);
       vis.setTransformColor(mdl_larm,trans,[0 0 1 i]); 
    end
    vis.update();
    b = false;
    while(~b)
        [I, b] = vis.getScreenshot();
    end
    writeVideo(vid,I);
end
close(vid);
%% Fade in actual markers

vid_file = 'fade_mark.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=0:0.005:1
    vis.addMarker('L',l_bar_end,[0.6 0.6 0 i]);
    vis.addMarker('R',r_bar_end,[0.6 0.6 0 i]);
    vis.update();
    b = false;
    while(~b)
        [I, b] = vis.getScreenshot();
    end
    writeVideo(vid,I);
end
close(vid);

%% Set up and run Multi Chain EKF 

%Reset model
mdl_rarm.position(:) = 0;mdl_rarm.velocity(:) = 0;mdl_rarm.acceleration(:) = 0;
mdl_larm.position(:) = 0;mdl_larm.velocity(:) = 0;mdl_larm.acceleration(:) = 0;
mdl_rarm.forwardPosition();
mdl_larm.forwardPosition();

ekf = MC_EKF_Q_DQ_DDQ();
ekf.addModel(mdl_larm);
ekf.addModel(mdl_rarm);
ekf.observation_noise = diag(repmat([0.1 0.1 0.1 1 1 1],1,6));
ekf.covariance = ekf.process_noise;
state_est_imu = zeros(size(dataInstance_imu_L_post.data,1),ekf.sizeX);

%Add Constraints between the two chains
% C1 = EKF_Constraint_T0EE('C1',mdl_larm,'frame_lhand_end',eye(4),mdl_rarm,'frame_rhand_end',eye(4));
% C1.type(:) = false;C1.type = logical([1 0 0 1 0 1]');
% ekf.addConstraint(C1);

vis = rlVisualizer('vis',1920,1080);
vis.addModel(mdl_rarm);
vis.addModel(mdl_larm);
%Set Constraint Colors of constraints
const_color = [0 1 1 0]';
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
vis.setViewportPose(viewport_init);
%Add markers
vis.addMarker('L',l_bar_end,[0.6 0.6 0 1]);
vis.addMarker('R',r_bar_end,[0.6 0.6 0 1]);
%Add imus
for j=1:numel(mdl_rarm.sensors)
    vis.setSensorColor(mdl_rarm,mdl_rarm.sensors(j),[1 0 1 1]);
    t_indx = find(contains({mdl_rarm.transforms.name},[mdl_rarm.sensors(j).name 't']),1);
    trans = mdl_rarm.transforms(t_indx);
    vis.setTransformColor(mdl_rarm,trans,[0 0 1 1]);
end
for j=1:numel(mdl_larm.sensors)
    vis.setSensorColor(mdl_larm,mdl_larm.sensors(j),[1 0 1 1]);
    t_indx = find(contains({mdl_larm.transforms.name},[mdl_larm.sensors(j).name 't']),1);
    trans = mdl_larm.transforms(t_indx);
    vis.setTransformColor(mdl_larm,trans,[0 0 1 1]);
end
vis.setViewportPose(viewport_init)
vis.update();

%% Pivot to top view 
viewport_pivot = [0.3940    0.0597         0   -0.6000   88.2000    2.5232]';
x = [1 200];
xq = 1:200;
v = [viewport_init viewport_pivot];
viewport_poses = zeros(6,200);
for i=1:6
    viewport_poses(i,:) = interp1(x,v(i,:),xq);
end

vid_file = 'pivot_top.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=1:size(viewport_poses,2)
   vis.setViewportPose(viewport_poses(:,i));
   vis.update();
   b = false;
   while(~b)
       [I, b] = vis.getScreenshot();
   end
   writeVideo(vid,I);
end
close(vid);

%% Visualize Constraints
vid_file = 'fade_const.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=0:0.005:1
    vis.setSensorColor(C.models(1),C.sensors(1),[0 1 1 i]');
    vis.setSensorColor(C.models(2),C.sensors(2),[0 1 1 i]');
    vis.update();
    b = false;
    while(~b)
        [I, b] = vis.getScreenshot();
    end
    writeVideo(vid,I);
end
close(vid);

%% Pivot back to front view
vid_file = 'pivot_front.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=size(viewport_poses,2):-1:1
   vis.setViewportPose(viewport_poses(:,i));
   vis.update();
   b = false;
   while(~b)
       [I, b] = vis.getScreenshot();
   end
   writeVideo(vid,I);
end
close(vid);

%% Run EKF and see what happens
dt = dataInstance_imu_L_post.dt;

f_l = mdl_larm.getFrameByName('frame_lhand_end');
f_r = mdl_rarm.getFrameByName('frame_rhand_end');

%Transformation from right EE to left EE in left EE frame
Tl_r = zeros(4,4,size(dataInstance_imu_L_post.data,1));

vid_file = 'slide_no_const.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=1:size(dataInstance_imu_L_post.data,1)
    ekf.run_iteration(dt,[dataInstance_imu_L_post.data(i,:) dataInstance_imu_R_post.data(i,:)]');
    state_est_imu(i,:) = ekf.state;
    
    Tl_r(:,:,i) = SE3.fastInverse(f_l.t)*f_r.t;
    l_mark = (trc.data.LEFT_U(i,:)+trc.data.LEFT_L(i,:))/2;
    r_mark = (trc.data.RIGHT_U(i,:)+trc.data.RIGHT_L(i,:))/2;
    vis.addMarker('L',l_mark,[0.6 0.6 0 1]);
    vis.addMarker('R',r_mark,[0.6 0.6 0 1]);
    vis.update();
    b = false;
    while(~b)
        [I, b] = vis.getScreenshot();
    end
    writeVideo(vid,I);
end
close(vid);

%% No constraint pivot to top 
viewport_top = [0.3940    0.0597         0   -0.6000   88.2000    2.5232]';
x = [1 200];
xq = 1:200;
v = [viewport_init viewport_top];
viewport_poses = zeros(6,200);
for i=1:6
    viewport_poses(i,:) = interp1(x,v(i,:),xq);
end


vid_file = 'slide_no_const_pivot_top.avi';
vid = VideoWriter([vid_path vid_file]);
open(vid);

for i=1:size(viewport_poses,2)
   vis.setViewportPose(viewport_poses(:,i));
   vis.update();
   b = false;
   while(~b)
       [I, b] = vis.getScreenshot();
   end
   writeVideo(vid,I);
end
close(vid);


%% Lets look at plots

%Translation of right end effector in left end effector frame
tmp = permute(Tl_r,[1,3,2]);
t_l_r = tmp(1:3,:,4)';

%X axis of EE 2 in frame of EE 1
x_l_r = tmp(1:3,:,1)';


