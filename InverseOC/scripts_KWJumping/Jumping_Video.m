%Generate Pretty Video of Jumping

%% Load completed EKF
%This is reading TRC files and stuff
addpath('C:\aslab\projects\AutoRehabSystem\Common');
addpath('C:\aslab\projects\AutoRehabSystem\Common\Classes');
EKFCodePath = 'C:\aslab\projects\vjoukov\General_FKEKF\DynamicsModelMatlab\MatlabWrapper\';
addpath(genpath(EKFCodePath));

Subj = 'P10';
dataFolder = ['C:\aslab\data\Jumping_To_Target_Data_Collection_2018\Jumping_Data_All_TRC\' Subj '_filtered\'];
template_dataName = [Subj '_target_85_1_1_clean-' Subj '_template'];
dataName = [Subj '_target_85_1_1-' Subj '_template'];

load([dataFolder dataName '.mat']);


%% Create model and visualizer 
%Create Euler Model
initPosWindow = [40, 100];
[mdl_eul, ~, initPos] = createJumpModel(template_trc, 1, initPosWindow, 0.06, EKFCodePath);

mdl_eul.position = ekf_clean_state(1,1:numel(mdl_eul.joints));
mdl_eul.forwardPosition;
% %Shift the Lie group model by 1 meter in x and y to visualize better
% mdl_lie.transforms(1).t(1:2,4) = [1 1]';
% mdl_lie.forwardPosition();

vis = rlVisualizer('vis',1280,720);
vis.addModel(mdl_eul);
vis.setBackgroundColour([1 1 1]);

%Set to initial pose visualization
init_pose =[0.5367   -0.1597         0  146.5999   35.2000    4.6616]';
vis.setViewportPose(init_pose);
vis.update;

%% Fade in to explain what are markers what are sensors

%Make Sensors invisible 
for j=1:numel(mdl_eul.sensors)
   vis.setSensorColor(mdl_eul,mdl_eul.sensors(j),[1 0 1 0]);     
end

%Make Markers Invisible 
z_eul = pos_mes_eul_keep_const(1,:);
mes_obj.setMesArray(z_eul);
for m=1:numel(mes_obj)
    m_pos = mes_obj(m).mes(1:3);
    vis.addMarker(mdl_eul.sensors(m).name,m_pos,[0.5 0.5 0 0]);
end

pause(2);
v = VideoWriter('marker_fade.avi');
open(v);
for i=1:50
    for j=1:numel(mdl_eul.sensors)
        vis.setSensorColor(mdl_eul,mdl_eul.sensors(j),[1 0 1 i/50]);
    end
    vis.update();
    writeVideo(v,vis.getScreenshot);
end

for i=1:50
    z_eul = pos_mes_eul_keep_const(1,:);
    mes_obj.setMesArray(z_eul);
    for m=1:numel(mes_obj)
        m_pos = mes_obj(m).mes(1:3);
        vis.addMarker(mdl_eul.sensors(m).name,m_pos,[0.5 0.5 0 i/50]);
    end
    vis.update();
    writeVideo(v,vis.getScreenshot);
end

close(v);


%% Make the video 


v = VideoWriter('jumping_video.avi');
open(v);


for i=1:size(ekf_clean_state,1)
   mdl_eul.position =  ekf_raw_state(i,1:numel(mdl_eul.joints));
   mdl_eul.forwardPosition();
   
   z_eul = pos_mes_eul_keep_const(i,:);
   mes_obj.setMesArray(z_eul);
   for m=1:numel(mes_obj)
       m_pos = mes_obj(m).mes(1:3);
       vis.addMarker(mdl_eul.sensors(m).name,m_pos,[0.5 0.5 0 1]);
   end
   vis.update();
   pause(dt);
end




