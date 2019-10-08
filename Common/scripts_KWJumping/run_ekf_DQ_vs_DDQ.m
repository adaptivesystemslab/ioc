%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
clear vis; close all;

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'20'};


markerMAE_all_ave_DQ = zeros(36, numel(partsToRun)); % MAE for all markers, averaged over all time, for selected jump
markerMAE_all_ave_DDQ = zeros(36, numel(partsToRun)); % MAE for all markers, averaged over all time, for selected jump
% markerMAE_marker_ave = zeros(2, 30, numel(i_part)); % MAE for each marker, averaged over all time, for selected jump


for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};


    data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));
    
    
        % Large Nested FOR Loops for each jump recording
        for i_targ = 1:3
            for i_set = 1:2
                for i_jump = 1:6
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
%                         [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_v2(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);
                        [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_setLinkLengths(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath,partNum);

                        mdl_eul.forwardPosition;


                        %% Run EKF
                        %Initialize model, set to initPos (approximate start position)
                        mdl_eul.position = initPos;
                        mdl_eul.velocity(:) = 0;
                        mdl_eul.acceleration(:) = 0;
                        mdl_eul.forwardPosition;
                        mdl_eul.forwardVelocity;


                        %Create the Measurement vector by concatenating marker data
                        mes_eul = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);

                        for i=1:numel(mdl_eul.sensors)
                            indxs = i*3-2:i*3;
                            mes_eul(:,indxs) = trcMod_eul.data.(mdl_eul.sensors(i).name);
                        end

                        mes_eul_predict = zeros(size(mes_eul));
                        
                        % copy model and parameters for second DDQ model
                        mdl_eul2 = mdl_eul;
                        mes_eul_predict2 = mes_eul_predict;
                        
                        ekf_eul = EKF_Q_DQ(mdl_eul);
                        ekf_eul2 = EKF_Q_DQ_DDQ(mdl_eul2);

                        ekf_eul.observation_noise = eye(ekf_eul.sizeZ) * 0.01;
                        ekf_eul.covariance = eye(ekf_eul.sizeX) * 1;
                        ekf_eul.process_noise = eye(ekf_eul.sizeX) * 1.05;

                        ekf_eul2.observation_noise = eye(ekf_eul2.sizeZ) * 0.01;
                        ekf_eul2.covariance = eye(ekf_eul2.sizeX) * 1;
                        ekf_eul2.process_noise = eye(ekf_eul2.sizeX) * 1.05;


                        % For new marker swapping MatlabWrapper
                        mes_obj = SensorMeasurement(1,30);
                        [mes_obj.size] = deal(3);
                        [mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);
                        matches = 1:numel(mdl_eul.sensors);
                        matches = [matches' matches'];


                        tic;
                        for i=1:size(mes_eul,1)


                            z_eul = mes_eul(i,:);
                            mes_obj.setMesArray(z_eul); % For marker swapping MatlabWrapper
                            u = dt;

                            ekf_eul.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper
                            mes_eul_predict(i,:) = ekf_eul.makeMeasure(ekf_eul.state).getMesArray(); % For marker swapping MatlabWrapper

                            ekf_eul2.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper
                            mes_eul_predict2(i,:) = ekf_eul2.makeMeasure(ekf_eul2.state).getMesArray(); % For marker swapping MatlabWrapper

                        end
                        EKF_runtime = toc

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
                            mae_eul_FB(i) = vpa(mean(norm_err_eul(i,:)));
                        end

                        % MAE for each marker, averaged over all time
                        mae_eul_markers = zeros(size(norm_err_eul,2),1);
                        for i = 1:size(norm_err_eul,2)
                            mae_eul_markers(i) = mean(norm_err_eul(:,i));
                        end
                        
                        
                        % Copy MAE calculations for second DDQ model
                        err_eul2 = mes_eul - mes_eul_predict2;

                        %Distance between each marker real and prediction, at each timestep
                        norm_err_eul2 = zeros(size(err_eul2,1),size(err_eul2,2)/3);
                        for i=1:(size(err_eul2,2)/3)
                            norm_err_eul2(:,i) = sqrt(sum(err_eul2(:,i*3-2:i*3).^2,2));
                        end

                        % MAE for full body (FB), at each timestep
                        mae_eul_FB2 = zeros(size(norm_err_eul2,1),1);
                        for i = 1:size(norm_err_eul2,1)
                            mae_eul_FB2(i) = vpa(mean(norm_err_eul2(i,:)));
                        end

                        % MAE for each marker, averaged over all time
                        mae_eul_markers2 = zeros(size(norm_err_eul2,2),1);
                        for i = 1:size(norm_err_eul2,2)
                            mae_eul_markers2(i) = mean(norm_err_eul2(:,i));
                        end
                        


                        currJump = 12*(i_targ-1) + 6*(i_set-1) + i_jump;
                        markerMAE_all_ave_DQ(currJump,i_part) = vpa(mean(mae_eul_FB));
                        markerMAE_all_ave_DDQ(currJump,i_part) = vpa(mean(mae_eul_FB2));
%                         markerMAE_marker_ave(EKF_procNoise,EKF_obsNoise,i_part,:) = mae_eul_markers;
                        
%                         figure(1); 
%                         for i_fig = 1:30
%                             clf; hold on; grid on;
%                             plot(norm_err_eul(:,i_fig),'r');
%                             plot(norm_err_eul2(:,i_fig),'b');
%                             legend('DQ','DDQ');
%                             title(['marker num ' num2str(i_fig)]);
%                             w = waitforbuttonpress();
%                         end
                        
                        
                        disp([dataName ' filtering complete'])
                    end
        % END MAIN NESTED FOR LOOPS
                end
            end
        end
                
                
                
        % End of EKF parameter optimization nested for loops

    % finished this participant
end

