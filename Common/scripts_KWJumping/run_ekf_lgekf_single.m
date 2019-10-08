%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
<<<<<<< .mine
clear vis; close all;

||||||| .r5286
clear all; close all;

=======
>>>>>>> .r5364
addpath('..\..\');
<<<<<<< .mine
% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
||||||| .r5286
% addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking'));

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
=======

%This is reading TRC files and stuff
addpath('C:\aslab\projects\AutoRehabSystem\Common');
addpath('C:\aslab\projects\AutoRehabSystem\Common\Classes');

%ADD MAGIC EKF STUFF
EKFCodePath = 'C:\aslab\projects\vjoukov\General_FKEKF\DynamicsModelMatlab\MatlabWrapper\';
>>>>>>> .r5364
addpath(genpath(EKFCodePath));

Subj = 'P10';
dataFolder = ['C:\aslab\data\Jumping_To_Target_Data_Collection_2018\Jumping_Data_All_TRC\' Subj '_filtered\'];
template_dataName = [Subj '_target_85_1_1_clean-' Subj '_template'];
dataName = [Subj '_target_85_1_1-' Subj '_template'];

<<<<<<< .mine
dataFolder = 'Vlad_EKF_test/';
template_dataName = 'P01_target_85_2_6_clean-P01_template';
dataName = 'P01_target_85_2_6-P01_template';

partNum = '03';
[age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
S2sc = S2sc/100;

||||||| .r5364
dataFolder = 'Vlad_EKF_test/';
template_dataName = 'P01_target_85_2_6_clean-P01_template';
dataName = 'P01_target_85_2_6-P01_template';

partNum = '03';
[gender,height,weight,S2sc] = getPartData(partNum);
S2sc = S2sc/100;

=======
>>>>>>> .r5412
visualize = 1;

%% Read in marker data, make models, visualize
template_trc = readTrc([dataFolder template_dataName '.trc']);
trc = readTrc([dataFolder dataName '.trc']);
markerNames = fieldnames(template_trc.data);
markerNames = markerNames(3:end); % first 2 names are Frame # and Time
% Rotate markers so subject jumps in positive X direction
for m = 1:length(markerNames)
    template_trc.data.(markerNames{m}) = (rotz(-pi/2)*template_trc.data.(markerNames{m})')';
    template_trc.data.(markerNames{m})(template_trc.data.(markerNames{m})==0) = -1000;
    trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')'/1000;
end

<<<<<<< .mine
%Determine TO and landing frames in recording
footSpeedX = diff(trc.data.FRMed(:,1));
footTOMoveFrame = find(footSpeedX>=1,1);
footTOPosZ = mean(trc.data.FRMed(1:footTOMoveFrame,3));
footTOFrame = find(trc.data.FRMed(:,3)>(footTOPosZ+10),1);
startFrame = footTOFrame - 100; %add 1 second to start
landPos = (44*25.4)*0.55 - 4*25.4;
% landPos = (mean(trcSL.data.LAND_BL(1:200,1)) + mean(trcSL.data.LAND_BR(1:200,1)))/2; %x-coord of front edge of target platform
footOverTargetFrame = find(trc.data.FRMed(:,1) >= landPos,1);
footLandFrame = footOverTargetFrame + find(footSpeedX(footOverTargetFrame:end)<=1,1);
endFrame = footLandFrame + trc.DataRate; %add 1 second to end


%Create Euler and Lie Group Models
initPosWindow = [50, 150];
[mdl_eul, trc, initPos] = createJumpModel(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);
% [mdl_lie, trc, initPos] = createJumpModel_Lie(trc, 1, [1, 100], 0.04, 0.06);

||||||| .r5286
%Determine TO and landing frames in recording
footSpeedX = diff(trc.data.FRMed(:,1));
footTOMoveFrame = find(footSpeedX>=1,1);
footTOPosZ = mean(trc.data.FRMed(1:footTOMoveFrame,3));
footTOFrame = find(trc.data.FRMed(:,3)>(footTOPosZ+10),1);
startFrame = footTOFrame - 100; %add 1 second to start
landPos = (44*25.4)*0.55 - 4*25.4;
% landPos = (mean(trcSL.data.LAND_BL(1:200,1)) + mean(trcSL.data.LAND_BR(1:200,1)))/2; %x-coord of front edge of target platform
footOverTargetFrame = find(trc.data.FRMed(:,1) >= landPos,1);
footLandFrame = footOverTargetFrame + find(footSpeedX(footOverTargetFrame:end)<=1,1);
endFrame = footLandFrame + trc.DataRate; %add 1 second to end


%Create Euler and Lie Group Models
initPosWindow = [50, 150];
[mdl_eul, trc, initPos] = createJumpModel(trc, 1, initPosWindow, S2sc, EKFCodePath);
% [mdl_lie, trc, initPos] = createJumpModel_Lie(trc, 1, [1, 100], 0.04, 0.06);

=======
%Create Euler Model
initPosWindow = [40, 100];
[mdl_eul, ~, initPos] = createJumpModel(template_trc, 1, initPosWindow, 0.06, EKFCodePath);
>>>>>>> .r5364
mdl_eul.forwardPosition;
% %Shift the Lie group model by 1 meter in x and y to visualize better
% mdl_lie.transforms(1).t(1:2,4) = [1 1]';
% mdl_lie.forwardPosition();

if(visualize)
    vis = rlVisualizer('vis', 640, 960);
    vis.addModel(mdl_eul);
    vis.update;
end;

%% Setup EKF
%Initialize model, set to initPos (approximate start position)
mdl_eul.position(:) = 0;

joint_names = {mdl_eul.joints.name};
mdl_eul.joints(find(ismember(joint_names,'lshoulder_jAbduction')==1)).position = pi/2
mdl_eul.joints(find(ismember(joint_names,'rshoulder_jAbduction')==1)).position = pi/2
mdl_eul.joints(find(ismember(joint_names,'lshoulder_jExtRotation')==1)).position = pi/2
mdl_eul.joints(find(ismember(joint_names,'rshoulder_jExtRotation')==1)).position = pi/2

%mdl_eul.joints(find(ismember(joint_names,'relbow_jFlexion')==1)).position = pi/2
%mdl_eul.joints(find(ismember(joint_names,'lelbow_jFlexion')==1)).position = pi/2

mdl_eul.velocity(:) = 0;
mdl_eul.acceleration(:) = 0;
mdl_eul.forwardPosition;
mdl_eul.forwardVelocity;

if(visualize)
    vis.update;
end

%Create the Measurement vector by concatenating marker data
pos_mes_eul = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);
pos_mes_eul_template = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);
vel_mes_eul = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);
mes_eul = [pos_mes_eul];
for i=1:numel(mdl_eul.sensors)
    indxs = i*3-2:i*3;
    pos_mes_eul(:,indxs) = trc.data.(mdl_eul.sensors(i).name);
    pos_mes_eul_template(:,indxs) = template_trc.data.(mdl_eul.sensors(i).name)/1000;
    vel_mes_eul(:,indxs) = [reshape(smooth(diff(trc.data.(mdl_eul.sensors(i).name))),[],3)*trc.DataRate;[0 0 0]];
    %mes_eul(:,i*6-5:i*6) = [pos_mes_eul(:,indxs) vel_mes_eul(:,indxs)];
    pos_mes_this_marker = trc.data.(mdl_eul.sensors(i).name);
    missing_markers = pos_mes_this_marker==0;
    pos_mes_this_marker(missing_markers) = 0;
    pos_mes_eul(:,indxs) = pos_mes_this_marker;
    mes_eul(:,i*3-2:i*3) = pos_mes_eul(:,indxs);
end

%Create a measurement that keeps position constant 
pos_mes_eul_keep_const = pos_mes_eul;
for i=1:numel(mdl_eul.sensors)
    indxs = i*3-2:i*3;
    pos_mes_eul(:,indxs) = trc.data.(mdl_eul.sensors(i).name);
    pos_mes_this_marker = trc.data.(mdl_eul.sensors(i).name);
    missing_markers = ismember(pos_mes_this_marker,[0 0 0],'rows');
    %Note this actually finds the index right before missing
    missing_starts = strfind(missing_markers',[0 1]);
    missing_ends = strfind(missing_markers',[1 0]);
    if(numel(missing_starts) > numel(missing_ends))
        missing_ends = [missing_ends size(pos_mes_this_marker,1)];
    end
    if(missing_markers(1))
        missing_starts = [1 missing_starts];
    end
    for(k=1:numel(missing_starts))
        pos_mes_this_marker(missing_starts(k):missing_ends(k),:) = ...
            repmat(pos_mes_this_marker(missing_starts(k),:),missing_ends(k)-missing_starts(k)+1,1);
    end
    pos_mes_eul_keep_const(:,indxs) = pos_mes_this_marker;
end

%Time Difference between samples
dt = 1/trc.DataRate;

mes_eul_predict = zeros(size(mes_eul));

ekf_eul = GP_EKF_Q_DQ(mdl_eul);

%EKF that will use clean data
ekf_clean = GP_EKF_Q_DQ(mdl_eul);

%EKF that will use raw data but set marker noise to super high when marker
%is missing
ekf_raw = GP_EKF_Q_DQ(mdl_eul);

%Event Handling for markers
event_handler = EventHandler();
ekf_eul.addlistener('matched_callback',@event_handler.matched_marker);
ekf_eul.addlistener('missing_callback',@event_handler.missing_marker);
ekf_eul.addlistener('swap_callback',@event_handler.swap_marker);

%Put constraints on the EKF
%Elbow Constraints
%ekf_eul.lower_state_bound(find(ismember(joint_names,'lshoulder_jExtRotation')==1)) = -pi/4;
%ekf_eul.upper_state_bound(find(ismember(joint_names,'lshoulder_jElevation')==1)) = pi;

%FOR POSITION ONLY
ekf_eul.observation_noise = diag(repmat([0.01 0.01 0.01],1,numel(mdl_eul.sensors)));
ekf_clean.observation_noise = ekf_eul.observation_noise;
ekf_raw.observation_noise = ekf_eul.observation_noise;
ekf_eul.covariance = eye(ekf_eul.sizeX) * 1;
ekf_clean.covariance = eye(ekf_eul.sizeX) * 1;
ekf_raw.covariance = eye(ekf_eul.sizeX) * 1;

%Calculating Process Noise
eta = 1000;
dim = numel(mdl_eul.joints);
%Accel noise
%G = [ones(dim,1)*dt^2/2*eta; ones(dim,1)*dt*eta; ones(dim,1)*eta];
%Velocity noise
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
ekf_eul.process_noise = P;
ekf_clean.process_noise = P;
ekf_raw.process_noise = P;

%Add some more to position so we dont do weird snapies
%ekf_eul.process_noise(1:dim,1:dim) = ekf_eul.process_noise(1:dim,1:dim) + 0.1;
%ekf_eul.process_noise(dim+1:dim*2,dim+1:dim*2) = ekf_eul.process_noise(dim+1:dim*2,dim+1:dim*2) + 0.1;

<<<<<<< .mine
% For new marker swapping MatlabWrapper
mes_obj = SensorMeasurement(1,30);
[mes_obj.size] = deal(3);
[mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);
matches = 1:numel(mdl_eul.sensors);
matches = [matches' matches'];

||||||| .r5286
=======
%This is the SensorMeasurement object array we need
mes_obj = SensorMeasurement(1,numel(markerNames));
%Set all of the measurement to have size 6 (position and velocity)
[mes_obj.size] = deal(3);
%Set all of the measurement to be markers
[mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);

ekf_state = zeros(size(mes_eul,1),numel(ekf_eul.state));


%Used to match arms in first few frames
sensor_names = {mdl_eul.sensors.name}';
%Define what the know matches are at the beginning
matchingSensors = {'ELLat','ELMed','WLLat','WLMed','ERLat','ERMed','WRMed','WRLat','WRMed','SR','SL','Neck','HL','HR'};

matches = zeros(numel(matchingSensors),2);
for i=1:numel(matchingSensors)
    matches(i,1) = find(ismember(sensor_names,matchingSensors{i})==1);
    matches(i,2) = find(ismember(sensor_names,matchingSensors{i})==1);
end


%Run first frame 150 times to converge
for i=1:100
    z_eul = pos_mes_eul_template(1,:);
    mes_obj.setMesArray(z_eul);
    ekf_eul.run_iteration(dt,mes_obj,matches);
      
    for m=1:numel(mes_obj)
        m_pos = mes_obj(m).mes(1:3);
        vis.addMarker(mdl_eul.sensors(m).name,m_pos);
    end
    vis.update();
end

init_model_pos = mdl_eul.position;

%% RUN Magic EKF

%Reset EKF initial covariances and stuff 
ekf_eul.covariance = ekf_clean.covariance;

%reset event handler
event_handler.mach_events = 0;
event_handler.swap_events = 0;
event_handler.miss_events = 0;

swap_frames = [];
miss_frames = [];
mach_frames = [];

>>>>>>> .r5364
tic;
for i=1:size(mes_eul,1)
<<<<<<< .mine
    
%     if(i==footTOFrame)
%         disp('footTOFrame, PRESS ANY KEY');
%         pause();
%     end
%     if(i==footLandFrame)
%         disp('footLandFrame, PRESS ANY KEY');
%         pause();
%     end
    
    z_eul = mes_eul(i,:)';
    mes_obj.setMesArray(z_eul); % For marker swapping MatlabWrapper
%     z_lie = mes_lie(i,:)';
    
||||||| .r5286
    
%     if(i==footTOFrame)
%         disp('footTOFrame, PRESS ANY KEY');
%         pause();
%     end
%     if(i==footLandFrame)
%         disp('footLandFrame, PRESS ANY KEY');
%         pause();
%     end
    
    z_eul = mes_eul(i,:)';
%     z_lie = mes_lie(i,:)';
    
=======
    z_eul = mes_eul(i,:);
    mes_obj.setMesArray(z_eul);
>>>>>>> .r5364
    u = dt;
    
<<<<<<< .mine
<<<<<<< .mine
    ekf_eul.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper
%     ekf_lie.run_iteration(u,z_lie);
||||||| .r5286
    ekf_eul.run_iteration(u,z_eul);
%     ekf_lie.run_iteration(u,z_lie);
=======
    %In the first few frames we force the matches to snap gooder
    if(i<80)
        ekf_eul.run_iteration(u,mes_obj,matches);
    else
        ekf_eul.run_iteration(u,mes_obj);
||||||| .r5364
    %In the first few frames we force the matches to snap gooder
    if(i<80)
        ekf_eul.run_iteration(u,mes_obj,matches);
    else
        ekf_eul.run_iteration(u,mes_obj);
=======
    if i == 1070
       disp('bp'); 
>>>>>>> .r5412
    end
    
    pre_swap = event_handler.swap_events;
    pre_mis = event_handler.miss_events;
    pre_match = event_handler.mach_events;
    ekf_eul.run_iteration(u,mes_obj);
    if(pre_swap ~= event_handler.swap_events)
       swap_frames(end+1) = i; 
    end
    if(pre_mis ~= event_handler.miss_events)
       miss_frames(end+1) = i; 
    end
    if(pre_match ~= event_handler.mach_events)
       mach_frames(end+1) = i; 
    end
    
    
    mes_eul_predict(i,:) = ekf_eul.makeMeasure(ekf_eul.state).getMesArray;
<<<<<<< .mine
>>>>>>> .r5364
    
||||||| .r5364
    
=======
>>>>>>> .r5412
    ekf_state(i,:) = ekf_eul.state;
    
<<<<<<< .mine
<<<<<<< .mine
    mes_eul_predict(i,:) = ekf_eul.makeMeasure(ekf_eul.state).getMesArray(); % For marker swapping MatlabWrapper
%     mes_lie_predict(i,:) = ekf_lie.makeMeasure(ekf_lie.state);

    jointAngles_eul = [jointAngles_eul; mdl_eul.position'];
%     jointAngles_lie = [jointAngles_lie; mdl_lie.position'];
    
%     mes_eul_predict_V(i,:) = vertcat(ekf_eul.model_handle.sensors.velocity)';
%     mes_lie_predict_V(i,:) = vertcat(ekf_lie.model_handle.sensors.velocity)';
%     
%     % Rotate marker angular and linear velocities back to world frame
%     for j=1:numel(ekf_eul.model_handle.sensors)
%         sens_eul =  ekf_eul.model_handle.sensors(j);
%         sens_lie =  ekf_lie.model_handle.sensors(j);
%         
%         Rlie = sens_lie.transform(1:3,1:3);
%         V1 = Rlie*sens_lie.velocity(1:3);
%         V2 = Rlie*sens_lie.velocity(4:6);
%         mes_lie_predict_V_rot(i,j*6-5:j*6) = [V1;V2];
%         
%         Reul = sens_eul.transform(1:3,1:3);
%         V1 = Reul*sens_eul.velocity(1:3);
%         V2 = Reul*sens_eul.velocity(4:6);
%         mes_eul_predict_V_rot(i,j*6-5:j*6) = [V1;V2];
%         
%     end
    
%     ekf_covariance(:,:,i) = ekf_eul.covariance;
%     lg_ekf_covariance(:,:,i) = ekf_lie.covariance;
    
    if(visualize)
||||||| .r5286
    mes_eul_predict(i,:) = ekf_eul.makeMeasure(ekf_eul.state);
%     mes_lie_predict(i,:) = ekf_lie.makeMeasure(ekf_lie.state);

    jointAngles_eul = [jointAngles_eul; mdl_eul.position'];
%     jointAngles_lie = [jointAngles_lie; mdl_lie.position'];
    
%     mes_eul_predict_V(i,:) = vertcat(ekf_eul.model_handle.sensors.velocity)';
%     mes_lie_predict_V(i,:) = vertcat(ekf_lie.model_handle.sensors.velocity)';
%     
%     % Rotate marker angular and linear velocities back to world frame
%     for j=1:numel(ekf_eul.model_handle.sensors)
%         sens_eul =  ekf_eul.model_handle.sensors(j);
%         sens_lie =  ekf_lie.model_handle.sensors(j);
%         
%         Rlie = sens_lie.transform(1:3,1:3);
%         V1 = Rlie*sens_lie.velocity(1:3);
%         V2 = Rlie*sens_lie.velocity(4:6);
%         mes_lie_predict_V_rot(i,j*6-5:j*6) = [V1;V2];
%         
%         Reul = sens_eul.transform(1:3,1:3);
%         V1 = Reul*sens_eul.velocity(1:3);
%         V2 = Reul*sens_eul.velocity(4:6);
%         mes_eul_predict_V_rot(i,j*6-5:j*6) = [V1;V2];
%         
%     end
    
%     ekf_covariance(:,:,i) = ekf_eul.covariance;
%     lg_ekf_covariance(:,:,i) = ekf_lie.covariance;
    
    if(visualize)
=======
    if(visualize && mod(i,3) == 0)
>>>>>>> .r5364
        for m=1:numel(z_eul)/3
            m_pos = z_eul(m*3-2:m*3);
||||||| .r5364
    if(visualize && mod(i,3) == 0)
        for m=1:numel(z_eul)/3
            m_pos = z_eul(m*3-2:m*3);
=======
    if(visualize && mod(i,100) == 0)
        for m=1:numel(mes_obj)
            m_pos = mes_obj(m).mes(1:3);
>>>>>>> .r5412
            vis.addMarker(mdl_eul.sensors(m).name,m_pos);
        end
        
        vis.update();
    end
    %disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
end
toc;

%% Now run the RAW DATA EKF that sets marker noise super high when markers are missing

%Reset EKF
mdl_eul.position = init_model_pos;
mdl_eul.velocity(:) = 0;
mdl_eul.forwardPosition();
mdl_eul.forwardVelocity();
ekf_raw.state(1:end/2) = init_model_pos;
ekf_raw.state(end/2+1:end) = 0;



%For clean EKF we force all 1-1 matches
matches = 1:numel(mdl_eul.sensors);
matches =[matches; matches]';
ekf_raw_state = zeros(size(mes_eul,1),numel(ekf_raw.state));
mes_raw_predict = zeros(size(mes_eul));

for i=1:size(mes_eul,1)
    z_eul = pos_mes_eul_keep_const(i,:);
    mes_obj.setMesArray(z_eul);
    u = dt;
    %For clean EKF we force all 1-1 matches
    ekf_raw.run_iteration(u,mes_obj,matches);
    mes_raw_predict(i,:) = ekf_raw.makeMeasure(ekf_raw.state).getMesArray;
    ekf_raw_state(i,:) = ekf_raw.state;
    if(visualize && mod(i,100) == 0)
        for m=1:numel(mes_obj)
            m_pos = mes_obj(m).mes(1:3);
            vis.addMarker(mdl_eul.sensors(m).name,m_pos);
        end
        vis.update();
    end
    disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
end


%% Now run the CLEAN DATA EKF

mdl_eul.position = init_model_pos;
mdl_eul.velocity(:) = 0;
mdl_eul.forwardPosition();
mdl_eul.forwardVelocity();
ekf_clean.state(1:end/2) = init_model_pos;

%For clean EKF we force all 1-1 matches
matches = 1:numel(mdl_eul.sensors);
matches =[matches; matches]';
ekf_clean_state = zeros(size(mes_eul,1),numel(ekf_eul.state));
mes_clean_predict = zeros(size(mes_eul));

for i=1:size(mes_eul,1)
    z_eul = pos_mes_eul_template(i,:);
    mes_obj.setMesArray(z_eul);
    u = dt;
    %For clean EKF we force all 1-1 matches
    ekf_clean.run_iteration(u,mes_obj,matches);
    mes_clean_predict(i,:) = ekf_clean.makeMeasure(ekf_clean.state).getMesArray;
    ekf_clean_state(i,:) = ekf_clean.state;
    if(visualize && mod(i,100) == 0)
        for m=1:numel(mes_obj)
            m_pos = mes_obj(m).mes(1:3);
            vis.addMarker(mdl_eul.sensors(m).name,m_pos);
        end
        vis.update();
    end
    disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
end


%% Compute RMSE Between Cleaned and EKF

markerRmse = zeros(numel(markerNames),1);
markerRmseClean = zeros(numel(markerNames),1);
markerRmseRaw = zeros(numel(markerNames),1);

for i=1:numel(markerNames)
    markerRmse(i) = sqrt(mean(sum((mes_eul_predict(100:end,i*3-2:i*3) - pos_mes_eul_template(100:end,i*3-2:i*3)).^2,2)));
    markerRmseClean(i) = sqrt(mean(sum((mes_clean_predict(100:end,i*3-2:i*3) - pos_mes_eul_template(100:end,i*3-2:i*3)).^2,2)));
    markerRmseRaw(i) = sqrt(mean(sum((mes_raw_predict(100:end,i*3-2:i*3) - pos_mes_eul_template(100:end,i*3-2:i*3)).^2,2)));
end

%% Pretty example plot of marker matching

marker = 'WRLat';
sensor_indx = find(ismember(sensor_names,marker)==1);
mis_markers = trc.data.WRLat(:,1) == 0;
time = 0:dt:size(mes_eul,1)*dt-dt;

mis_starts = strfind(mis_markers',[0 1]);
if(mis_markers(1) == 1)
    mis_starts = [1 mis_starts];
end

mis_ends = strfind(mis_markers',[1 0]);
no_mis_starts = [1 mis_ends]+1;
no_mis_ends = [mis_starts numel(time)];

colors = distinguishable_colors(10);
clf;
hold on;
p_gt_clean = plot(time,pos_mes_eul_template(:,sensor_indx*3-2),'color',[colors(1,:) 0.5],'LineWidth',10)
for i=1:numel(no_mis_starts)
    p_gt_raw = plot(time(no_mis_starts(i):no_mis_ends(i)),...
        trc.data.WLLat(no_mis_starts(i):no_mis_ends(i),1),'color',colors(2,:),'LineWidth',4)
end
p_raw = plot(time,mes_eul_predict(:,sensor_indx*3-2),'--','color',colors(3,:),'LineWidth',3)
p_clean = plot(time,mes_clean_predict(:,sensor_indx*3-2),'--','color',colors(6,:),'LineWidth',3)
set(gca, 'FontSize', 22)
ylabel('Position (m)','FontSize', 24);
xlabel('Time (sec)', 'FontSize', 24);
leg = legend([p_gt_clean p_gt_raw p_raw p_clean],'GT_{Clean}','GT_{Raw}','EKF_{Raw}','EKF_{Clean}')




%% Clean up and save
clear ekf_eul
clear ekf_clean
clear ekf_raw
clear vis
clear mdl_eul

save([dataFolder dataName '.mat']);






