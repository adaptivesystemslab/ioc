%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
clear vis; close all;

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);

visualize = 0;
plotJointAngles = 0;
pauseActive = 0;
saveMarkerData = 1;
saveJointAngles = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
%               '13','14','15','16','17','18','19','20','21','22'};
partsToRun = {'02','04','11','19'};

processNoise = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.4, 1.8, 2.5, 3.5, 5.0, 7.0, 10];
observationNoise = [0.001, 0.003, 0.006, 0.01, 0.02, 0.03, 0.05, 0.1];

markerMAE_all_ave = zeros(numel(processNoise), numel(observationNoise), numel(i_part)); % MAE for all markers, averaged over all time, for selected jump
markerMAE_marker_ave = zeros(numel(processNoise), numel(observationNoise), numel(i_part), 30); % MAE for each marker, averaged over all time, for selected jump


for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};


    data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));
    
    
    for EKF_procNoise = 1:numel(processNoise)
        for EKF_obsNoise = 1:numel(observationNoise)
%             for EKF_cov = 1:numel(covariance)
                
    
    
                % Large Nested FOR Loops for each jump recording
                for i_targ = 2%1:3
                    for i_set = 2%1:2
                        for i_jump = 3%1:6
                %             dataName = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                %                         '_' num2str(i_jump) '_clean-P' partNum '_template'];
                            dataName = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                                        '_' num2str(i_jump) '_clean-P']; %because some files names "P03-template"
                            file_ind = strfind(data_files,dataName);
                            file_ind = find(not(cellfun('isempty', file_ind)));
                            dataName = data_files{file_ind};

                            dataNameSL = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                                        '_' num2str(i_jump) '_clean-start_land'];

                %             if(exist(['../data/' dataFolder dataName '.trc'],'file')~=2)
                %                 disp([dataName '.TRC does not exist'])
                            if(exist(['../data/' dataFolder dataName],'file')~=2)
                                disp([dataName '.TRC does not exist'])

                            else
                                %% Read in marker data, make models, visualize
                %                 trc = readTrc(['../data/' dataFolder dataName '.trc']);
                                trc = readTrc(['../data/' dataFolder dataName]);
                                markerNames = fieldnames(trc.data);
                                markerNames = markerNames(3:end); % first 2 names are Frame # and Time
                                % Rotate markers so subject jumps in positive X direction
                                for m = 1:length(markerNames)
                                    trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
                                end
                %                 trc = fillInTRCFrames(trc);


                                % Read in start line and landing platform marker data
                                trcSL = readTrc(['../data/' dataFolder dataNameSL '.trc']);
                                markerNamesSL = fieldnames(trcSL.data);
                                markerNamesSL = markerNamesSL(3:end); % first 2 names are Frame # and Time
                                % Rotate markers so subject jumps in positive X direction
                                for m = 1:length(markerNamesSL)
                                    trcSL.data.(markerNamesSL{m}) = (rotz(-pi/2)*trcSL.data.(markerNamesSL{m})')';
                                end
                                trcSL = fillInTRCFrames(trcSL);

                                locationStart = double([mean(trcSL.data.START_L); mean(trcSL.data.START_R)])/1000; % convert to meters
                                locationLand = double([mean(trcSL.data.LAND_BL); mean(trcSL.data.LAND_BR);...
                                            mean(trcSL.data.LAND_TL); mean(trcSL.data.LAND_TR)])/1000;

                                %Determine TO and landing frames in recording
                                dt = 1/(trc.DataRate);                
                                footSpeedX = diff(trc.data.FRMed(:,1));
                                [~,midJumpFrame] = max(footSpeedX);
                                footTOFrame = find(footSpeedX(1:midJumpFrame)<=1,1,'last'); % equivalent to (1/1000)/dt = 0.2 [m/s]
                                footLandFrame = midJumpFrame + find(footSpeedX(midJumpFrame:end)<=1,1,'first');
                                startFrame = footTOFrame - trc.DataRate; %add 1 second to start
                                endFrame = footLandFrame + trc.DataRate; %add 1 second to end



                                %Create Euler and Lie Group Models
                                initPosWindow = [251,450];
                    %             initPosWindow = [(startFrame-100), startFrame];
            %                     [mdl_eul, trcMod_eul, initPos] = createJumpModel(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);
                                [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_v2(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);
                %                 [mdl_lie, trcMod_lie, initPos_lie] = createJumpModel_Lie(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);

                                mdl_eul.forwardPosition;
                                %Shift the Lie group model by 1 meter in x and y to visualize better
                %                 mdl_lie.transforms(1).t(1:2,4) = mdl_eul.transforms(1).t(1:2,4) + [1 1]';
                %                 mdl_lie.forwardPosition();

                                if(visualize)
                                    vis = rlVisualizer('vis', 640, 960);
                                    vis.addModel(mdl_eul);
                %                     vis.addModel(mdl_lie);
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

                %                 mdl_lie.position = initPos_lie;
                %                 mdl_lie.velocity(:) = 0;
                %                 mdl_lie.acceleration(:) = 0;
                %                 mdl_lie.forwardPosition;
                %                 mdl_lie.forwardVelocity;

                                if(visualize)
                                    vis.update;
                                end

                                %Create the Measurement vector by concatenating marker data
                                mes_eul = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);

                                for i=1:numel(mdl_eul.sensors)
                                    indxs = i*3-2:i*3;
                                    mes_eul(:,indxs) = trcMod_eul.data.(mdl_eul.sensors(i).name);
                                end
                                %Because of the shift in model placement need to add
                %                 mes_lie = mes_eul + repmat([1 1 0],size(mes_eul,1),numel(mdl_lie.sensors));


                                %Resample at a lower rate
                    %             resample = 1;
                                %Time Difference between samples
                    %             dt = 1/((trc.DataRate)/resample);

                    %             mes_eul = mes_eul(1:resample:end,:);
                                % mes_lie = mes_lie(1:resample:end,:);

                                mes_eul_predict = zeros(size(mes_eul));
                %                 mes_lie_predict = zeros(size(mes_lie));

                                % mes_eul_predict_V = zeros(size(mes_eul,1),2*size(mes_eul,2)); % marker velocity from ekf is 6x1 vector (not 3D position)
                                % mes_lie_predict_V = zeros(size(mes_lie,1),2*size(mes_lie,2));
                                % mes_eul_predict_V_rot = zeros(size(mes_eul,1),2*size(mes_eul,2)); % marker velocity from ekf is 6x1 vector (not 3D position)
                                % mes_lie_predict_V_rot = zeros(size(mes_lie,1),2*size(mes_lie,2));

%                                 ekf_eul = EKF_Q_DQ_DDQ(mdl_eul);
                                ekf_eul = EKF_Q_DQ(mdl_eul);
                %                 ekf_lie = LG_EKF_Q_DQ_DDQ(mdl_lie);

                                ekf_eul.observation_noise = eye(ekf_eul.sizeZ) * observationNoise(EKF_obsNoise);
                %                 ekf_lie.observation_noise = ekf_eul.observation_noise;

                                ekf_eul.covariance = eye(ekf_eul.sizeX) * 1;
                %                 ekf_lie.covariance = ekf_eul.covariance;

                                %Calculating Process Noise
                %                 eta = 300;
                %                 dim = ekf_lie.dim/3; % num markers
                %                 G = [ones(dim,1)*dt^2/2*eta; ones(dim,1)*dt*eta; ones(dim,1)*eta];
                %                 P_tmp = G*G';
                %                 P = zeros(size(P_tmp));
                %                 for i=1:3
                %                     for j=i:3
                %                         P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
                %                             diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim))) ;
                %                     end
                %                 end
                %                 P = P+P' - diag(diag(P));
                %                 ekf_eul.process_noise = P;
                %                 ekf_lie.process_noise = P;
                                ekf_eul.process_noise = eye(ekf_eul.sizeX) * processNoise(EKF_procNoise);
                %                 ekf_lie.process_noise = ekf_eul.process_noise;

                                ekf_covariance = zeros(size(ekf_eul.covariance,1),size(ekf_eul.covariance,2),size(mes_eul,1));
                %                 lg_ekf_covariance = zeros(size(ekf_lie.covariance,1),size(ekf_lie.covariance,2),size(mes_eul,1));

                                %We will also save hatInv of state instead of full out state
                                ekf_state = zeros(size(mes_eul,1),ekf_eul.sizeX);
                %                 lg_ekf_hatInvState = zeros(size(mes_eul,1),size(ekf_lie.hatinvState(ekf_lie.logState(ekf_lie.state)),1));

                                jointAngles_names = {mdl_eul.joints.name};
                                jointAngles_eul = [];
                %                 jointAngles_lie = [];


                                % For new marker swapping MatlabWrapper
                                mes_obj = SensorMeasurement(1,30);
                                [mes_obj.size] = deal(3);
                                [mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);
                                matches = 1:numel(mdl_eul.sensors);
                                matches = [matches' matches'];


                                tic;
                                for i=1:size(mes_eul,1)
                    %             for i=(startFrame-100):size(mes_eul,1)

                                %     if(i==footTOFrame)
                                %         disp('footTOFrame, PRESS ANY KEY');
                                %         pause();
                                %     end
                                %     if(i==footLandFrame)
                                %         disp('footLandFrame, PRESS ANY KEY');
                                %         pause();
                                %     end

                                    z_eul = mes_eul(i,:);
                                    mes_obj.setMesArray(z_eul); % For marker swapping MatlabWrapper
                %                     z_lie = mes_lie(i,:)';

                                    u = dt;

                                    ekf_eul.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper
                %                     ekf_lie.run_iteration(u,z_lie);

                                    ekf_state(i,:) = ekf_eul.state;
                %                     lg_ekf_hatInvState(i,:) = ekf_lie.hatinvState(ekf_lie.logState(ekf_lie.state));

                                    mes_eul_predict(i,:) = ekf_eul.makeMeasure(ekf_eul.state).getMesArray(); % For marker swapping MatlabWrapper
                %                     mes_lie_predict(i,:) = ekf_lie.makeMeasure(ekf_lie.state);

                                    jointAngles_eul = [jointAngles_eul; mdl_eul.position'];
                %                     jointAngles_lie = [jointAngles_lie; mdl_lie.position']; % NOTE: this doesn't work, since "mdl_lie.position" is in LieMapping
                %                     jointAngles_lie = getLieModelJAs(mdl_lie);


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

                                    ekf_covariance(:,:,i) = ekf_eul.covariance;
                %                     lg_ekf_covariance(:,:,i) = ekf_lie.covariance;


                                    if(visualize && (mod(i,5)==0) )
                %                         for m=1:numel(z_eul)/3
                %                             m_pos = z_eul(m*3-2:m*3);
                %                             vis.addMarker([num2str(m)],m_pos);
                %                         end

                %                         for m=1:numel(z_lie)/3
                %                             m_pos = z_lie(m*3-2:m*3);
                %                             vis.addMarker([num2str(m) 'L'],m_pos);
                %                         end

                                        vis.update();
                                        %pause(dt);
                                %         pause(0.03);
                                    end




                    %                 disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);

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
                                EKF_runtime = toc

                                if(visualize)
                                    clear vis;
                                end

                                %% Computing data and output

                                err_eul = mes_eul - mes_eul_predict;

                                %Distance between each marker real and prediction, at each timestep
                                norm_err_eul = zeros(size(err_eul,1),size(err_eul,2)/3);
                                for i=1:(size(err_eul,2)/3)
                                    norm_err_eul(:,i) = sqrt(sum(err_eul(:,i*3-2:i*3).^2,2));
                                end

                                % MAE for full body (FB), at each timestep
                                mae_eul_FB = zeros(size(norm_err_eul,1),1);
                                for i = 1:size(norm_err_eul,1)
                                    mae_eul_FB(i) = mean(norm_err_eul(i,:));
                                end

                                % MAE for each marker, averaged over all time
                                mae_eul_markers = zeros(size(norm_err_eul,2),1);
                                for i = 1:size(norm_err_eul,2)
                                    mae_eul_markers(i) = mean(norm_err_eul(:,i));
                                end
                                
                                
                                
                                markerMAE_all_ave(EKF_procNoise,EKF_obsNoise,i_part) = mean(mae_eul_FB);
                                markerMAE_marker_ave(EKF_procNoise,EKF_obsNoise,i_part,:) = mae_eul_markers;
                                
                                
                                
                                %% Save mes_### and mes_###_predict data

%                                 if(saveJointAngles)
%                                     mydir = pwd;
%                                     idcs = strfind(mydir,'\');
%                                     newdir = mydir(1:idcs(end)); 
%                                     saveFilePath = [newdir 'results\RESULTS_JA_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
% 
%                                     world2base = mdl_eul.transforms(1).t;
% 
%                                     save(saveFilePath,'jointAngles_names','jointAngles_eul',...
%                                             'footTOFrame','footLandFrame','startFrame','endFrame',...
%                                         'modelLinks','world2base','locationStart','locationLand');
%                 %                     save(saveFilePath,'jointAngles_names','jointAngles_eul','jointAngles_lie',...
%                 %                             'footTOFrame','footLandFrame','startFrame','endFrame');
% 
%                                     disp([dataName ' JA data saved'])
%                                 end


                                disp([dataName ' filtering complete'])
                            end
                % END MAIN NESTED FOR LOOPS
                        end
                    end
                end
                
                
                
                % End of EKF parameter optimization nested for loops
%             end
        end
    end

    % finished this participant
end

