%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
clear vis; close all;

addpath('..\..\');
% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
addpath(genpath(EKFCodePath));



dataFolder = 'Vlad_EKF_test/';
dataName = 'P03_target_55_1_1_clean-P03_template';

partNum = '03';
[age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
S2sc = S2sc/100;

visualize = 1;
pauseActive = 0;
plotJA = 0;

%% Read in marker data, make models, visualize
trc = readTrc(['../data/' dataFolder dataName '.trc']);
markerNames = fieldnames(trc.data);
markerNames = markerNames(3:end); % first 2 names are Frame # and Time
% Rotate markers so subject jumps in positive X direction
for m = 1:length(markerNames)
    trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
end

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

mdl_eul.forwardPosition;
% %Shift the Lie group model by 1 meter in x and y to visualize better
% mdl_lie.transforms(1).t(1:2,4) = [1 1]';
% mdl_lie.forwardPosition();

if(visualize)
    vis = rlVisualizer('vis', 640, 960);
    vis.addModel(mdl_eul);
    % vis.addModel(mdl_lie);
    vis.update;
end

% %Visualize the markers
% for m=1:numel(markerNames)
%     mp = trc_markers.data.(markerNames{m})(1,:);
%     %For some reason we have to call ADD marker twice before its placed
%     %correctly
%     vis.addMarker(markerNames{m},mp);
% %     vis.addMarker(m_names{m},mp);
% end
% vis.update();

if(pauseActive)
    disp('PRESS ANY KEY TO CONTINUE');
    pause();
end

%% Run EKF
%Initialize model, set to initPos (approximate start position)
mdl_eul.position = initPos;
mdl_eul.velocity(:) = 0;
mdl_eul.acceleration(:) = 0;
mdl_eul.forwardPosition;
mdl_eul.forwardVelocity;

% mdl_lie.position = initPos;
% mdl_lie.velocity(:) = 0;
% mdl_lie.acceleration(:) = 0;
% mdl_lie.forwardPosition;
% mdl_lie.forwardVelocity;

if(visualize)
    vis.update;
end

%Create the Measurement vector by concatenating marker data
mes_eul = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);

for i=1:numel(mdl_eul.sensors)
    indxs = i*3-2:i*3;
    mes_eul(:,indxs) = trc.data.(mdl_eul.sensors(i).name);
end
%Because of the shift in model placement need to add
% mes_lie = mes_eul + repmat([1 1 0],size(mes_eul,1),numel(mdl_lie.sensors));


%Resample at a lower rate
resample = 1;
%Time Difference between samples
dt = 1/(100/resample);

mes_eul = mes_eul(1:resample:end,:);
% mes_lie = mes_lie(1:resample:end,:);

mes_eul_predict = zeros(size(mes_eul));
% mes_lie_predict = zeros(size(mes_lie));

% mes_eul_predict_V = zeros(size(mes_eul,1),2*size(mes_eul,2)); % marker velocity from ekf is 6x1 vector (not 3D position)
% mes_lie_predict_V = zeros(size(mes_lie,1),2*size(mes_lie,2));
% mes_eul_predict_V_rot = zeros(size(mes_eul,1),2*size(mes_eul,2)); % marker velocity from ekf is 6x1 vector (not 3D position)
% mes_lie_predict_V_rot = zeros(size(mes_lie,1),2*size(mes_lie,2));

ekf_eul = EKF_Q_DQ_DDQ(mdl_eul);
% ekf_lie = LG_EKF_Q_DQ_DDQ(mdl_lie);

ekf_eul.observation_noise = eye(ekf_eul.sizeZ) * 0.01;
% ekf_lie.observation_noise = ekf_eul.observation_noise;

ekf_eul.covariance = eye(ekf_eul.sizeX) * 1;
% ekf_lie.covariance = ekf_eul.covariance;

%Calculating Process Noise
% eta = 300;
% dim = ekf_lie.dim/3; % num markers
% G = [ones(dim,1)*dt^2/2*eta; ones(dim,1)*dt*eta; ones(dim,1)*eta];
% P_tmp = G*G';
% P = zeros(size(P_tmp));
% for i=1:3
%     for j=i:3
%         P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
%             diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
%     end
% end
% P = P+P' - diag(diag(P));
% ekf_eul.process_noise = P;
% ekf_lie.process_noise = P;
ekf_eul.process_noise = eye(ekf_eul.sizeX) * 1;
% ekf_lie.process_noise = ekf_eul.process_noise;

% ekf_state_est = zeros(ekf_eul.sizeX,size(mes_eul,1));
% lg_ekf_state_est = zeros(ekf_lie.sizeX,size(mes_eul,1));

% ekf_covariance = zeros(size(ekf_eul.covariance,1),size(ekf_eul.covariance,2),size(mes_eul,1));
% lg_ekf_covariance = zeros(size(ekf_lie.covariance,1),size(ekf_lie.covariance,2),size(mes_eul,1));

%We will also save hatInv of state instead of full out state
% ekf_state = zeros(size(mes_eul,1),ekf_eul.sizeX);
% lg_ekf_hatInvState = zeros(size(mes_eul,1),size(ekf_lie.hatinvState(ekf_lie.logState(ekf_lie.state)),1));

jointAngles_names = {mdl_eul.joints.name};
jointAngles_eul = [];
jointAngles_lie = [];

% For new marker swapping MatlabWrapper
mes_obj = SensorMeasurement(1,30);
[mes_obj.size] = deal(3);
[mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);
matches = 1:numel(mdl_eul.sensors);
matches = [matches' matches'];

tic;
for i=1:size(mes_eul,1)
    
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
    
    u = dt;
    
    ekf_eul.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper
%     ekf_lie.run_iteration(u,z_lie);
    
%     ekf_state(i,:) = ekf_eul.state;
%     lg_ekf_hatInvState(i,:) = ekf_lie.hatinvState(ekf_lie.logState(ekf_lie.state));
    
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
        for m=1:numel(z_eul)/3
            m_pos = z_eul(m*3-2:m*3);
            vis.addMarker(['Meul' num2str(m)],m_pos);
        end
%         for m=1:numel(z_lie)/3
%             m_pos = z_lie(m*3-2:m*3);
%             vis.addMarker(['Mlie' num2str(m)],m_pos);
%         end
        
        vis.update();
%         pause(dt);
%         pause(0.03);
    end
    disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
    
%     if mod(i,50)==0
%         disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
%         for m=1:numel(z_eul)/3
%             m_pos = z_eul(m*3-2:m*3);
%             vis.addMarker(['Meul' num2str(m)],m_pos,[1 1 0 1]);
%         end
% %         for m=1:numel(z_lie)/3
% %             m_pos = z_lie(m*3-2:m*3);
% %             vis.addMarker(['Mlie' num2str(m)],m_pos,[1 1 0 1]);
% %         end
%         
%         vis.update();
%     end
end
toc;

if(visualize)
    clear vis;
end

%% Computing data and output

%Plot joint angles
if(plotJA)
    % joints: (1-6) p0,p1,p2,r0,r1,r2, 
    % (7-9) backFB, backAxial, backLateral, 
    % (10-14) rshldrElev, rshldrAbd, rshldrExtRot, relbowFlex, relbowSup
    % (15-19) lshldrElev, lshldrAbd, lshldrExtRot, lelbowFlex, lelbowSup
    % (20-22) rhipFlex, rhipAbd, rhipExtRot, 
    % (23-26) rkneeExtend, rkneeExtRot, rankleDorsi, ranklePron
    % (27-29) lhipFlex, lhipAbd, lhipExtRot, 
    % (30-33) lkneeExtend, lkneeExtRot, lankleDorsi, lanklePron
    plotJoints(1).a = [7,8,9]; %back
    plotJoints(2).a = [10,15]; %shldr elevation
    plotJoints(3).a = [11,16]; %shldr abduction
    plotJoints(4).a = [12,17]; %shldr ext. rotation
    plotJoints(5).a = [13,18]; %elbow flexion
    plotJoints(6).a = [20,27]; %hip flexion
    % plotJoints(7).a = [21,28]; %hip abduction
    % plotJoints(8).a = [22,29]; %hip ext. rotation
    plotJoints(7).a = [23,30]; %knee flexion 
    % plotJoints(10).a = [24,31]; %knee ext. rotation
    plotJoints(8).a = [25,32]; %ankle dorsiflexion
    plotJoints(9).a = [26,33]; %ankle pronation
    plotNames = {'Back','Shoulder Elev.','Shoulder Abd.','Shoulder Ext. Rot.',...
        'Elbow Flex.','Hip Flex.','Knee Flex.','Ankle Dorsiflex.','Ankle Pronate'};

    dataLength = size(jointAngles_eul,1);
    figure(1);
    for i = 1:numel(plotJoints)
        clf; hold on;
        leg = [];
        for j = 1:numel(plotJoints(i).a)
            jj = plotJoints(i).a(j);
            plot(1:dataLength,rad2deg(jointAngles_eul(:,jj)),'LineWidth',2);
            leg = [leg; jointAngles_names(jj)];
        end
        plot([footTOFrame, footTOFrame],[-180,180],'g','LineWidth',2);
        plot([footLandFrame, footLandFrame],[-180,180],'r','LineWidth',2);
        plot([startFrame, startFrame],[-180,180],'k--');
        plot([endFrame, endFrame],[-180,180],'k--');
        title(plotNames(i));
        legend(cellstr(leg));
        xlabel('Frame Number');
        ylabel('Angle [deg]');
        ylim([-120,120]);
        grid on;
        w = waitforbuttonpress;
    end
end
% 
% %Compute marker position errors
% err_eul = mes_eul - mes_eul_predict;
% % err_lie = mes_lie - mes_lie_predict;
% 
% %Distance between each marker real and prediction, at each timestep
% norm_err_eul = zeros(size(err_eul,1),numel(mdl_eul.sensors));
% % norm_err_lie = zeros(size(err_lie,1),numel(mdl_lie.sensors));
% 
% for i=1:numel(mdl_eul.sensors)
%     norm_err_eul(:,i) = sqrt(sum(err_eul(:,i*3-2:i*3).^2,2));
% %     norm_err_lie(:,i) = sqrt(sum(err_lie(:,i*3-2:i*3).^2,2));
% end
% 
% 
% %Mean distance error for each marker, between startFrame and endFrame
% mean_err_eul = mean(norm_err_eul);
% [~,idx] = sort(mean_err_eul);
% mean_err_eul_sorted_names = markerNames(idx);
% 
% %Root mean squared error of ALL markers, at each timestep
% rmse_eul = zeros(size(norm_err_eul,1),1);
% % rmse_lie = zeros(size(norm_err_lie,1),1);
% 
% for i = 1:size(norm_err_eul,1)
%     rmse_eul(i) = sqrt( sum(norm_err_eul(i,:).^2)/size(norm_err_eul,2) );
% %     rmse_lie(i) = sqrt( sum(norm_err_lie(i,:).^2)/size(norm_err_lie,2) );
% end


% ekf_pos_cov_trace = zeros(size(ekf_covariance,3),1);
% lg_ekf_pos_cov_trace = zeros(size(ekf_covariance,3),1);
%
% lg_ekf_dim_pos = 0;
% for i=1:numel(lg_ekf.structure)
%     lg_ekf_dim_pos = lg_ekf_dim_pos + lg_ekf.structure(i).dim;
% end
%
% ekf_vel_cov_trace = zeros(size(ekf_covariance,3),1);
% lg_ekf_vel_cov_trace = zeros(size(ekf_covariance,3),1);
%
% ekf_accel_cov_trace = zeros(size(ekf_covariance,3),1);
% lg_ekf_accel_cov_trace = zeros(size(ekf_covariance,3),1);
%
% ekf_cov_trace = zeros(size(ekf_covariance,3),1);
% lg_ekf_cov_trace = zeros(size(ekf_covariance,3),1);
%
%
% for i=1:numel(ekf_pos_cov_trace)
%     %Ignore some joints that are not working
%     ekf_pos_cov_trace(i) = trace(ekf_covariance(1:36,1:36,i));
%     lg_ekf_pos_cov_trace(i) = trace(lg_ekf_covariance(1:36,1:36,i));
%
%     ekf_vel_cov_trace(i) = trace(ekf_covariance(5:8,5:8,i));
%     lg_ekf_vel_cov_trace(i) = trace(lg_ekf_covariance(5:8,5:8,i));
%
%     ekf_accel_cov_trace(i) = trace(ekf_covariance(9:12,9:12,i));
%     lg_ekf_accel_cov_trace(i) = trace(lg_ekf_covariance(9:12,9:12,i));
%
%     ekf_cov_trace(i) = trace(ekf_covariance(:,:,i));
%     lg_ekf_cov_trace(i) = trace(lg_ekf_covariance(:,:,i));
%
% end



%% Figures
% 
% %Marker RMS Error Total Average
% figure(2);clf;
% plot(rmse_eul);
% hold on
% % plot(rmse_lie);
% % plot(rmse_cmu);
% % legend('EKF');
% title(['RMSE Of All Markers']);
% xlabel('Frame');
% ylabel('Distance [m]');
% xlim([0 size(rmse_eul,1)]);
% 
% %Individual Marker RMS Error
% figure(3); clf; hold on;
% lineColours = {'k','r','b','k','k','k','k','k',... %torso
%     'r','r','r','r',... %left arm
%     'b','b','b','b',... %right arm
%     'g','g','g','g','g','g','g',... %left leg
%     'm','m','m','m','m','m','m'}; %right leg
% lineStyles = {'-','-','-','--','--','-.','-.','-.',... %torso
%     '--','--','-.','-.',... %left arm
%     '--','--','-.','-.',... %right arm
%     '-','-','--','--','-.','-.','-.',... %left leg
%     '-','-','--','--','-.','-.','-.'}; %right leg
% for i = 1:length(markerNames)
%     plot(norm_err_eul(:,i),'Color',lineColours{i},'LineStyle',lineStyles{i});
% end
% legend(markerNames);
% title(['RMSE Of individual markers']);
% xlabel('Frame');
% ylabel('Distance [m]');
% axis([0,size(norm_err_eul,1),0,max(max(norm_err_eul))]);




