function f = tuneEKFParams_allParts(x)

%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'}; % all participants took 53 hours, procNoise_optim = 1.0364
% partsToRun = {'02','05','08','11','14','17','20'}; % these 7 participants took 11 hours, procNoise_optim = 2.2442
% partsToRun = {'08'};

markerMAE_all_ave = zeros(numel(partsToRun),36); % MAE for all markers, averaged over all time, for all jumps
markerMAE_marker_ave = zeros(numel(partsToRun),36,30); % MAE for each marker, averaged over all time, for all jumps

% tic;
for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
    [~,~,~,~,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};


    data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));
    
    % Large Nested FOR Loops for each jump recording
    for i_targ = 1:3
        for i_set = 1:2
            for i_jump = 1:6
%                 dataName = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
%                             '_' num2str(i_jump) '_clean-P' partNum '_template'];
                dataName = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                            '_' num2str(i_jump) '_clean-P']; %because some files names "P03-template"
                file_ind = strfind(data_files,dataName);
                file_ind = find(not(cellfun('isempty', file_ind)));
                dataName = data_files{file_ind};

%                 dataNameSL = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
%                             '_' num2str(i_jump) '_clean-start_land'];

%                 if(exist(['../data/' dataFolder dataName '.trc'],'file')~=2)
%                     disp([dataName '.TRC does not exist'])
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


%                     % Read in start line and landing platform marker data
%                     trcSL = readTrc(['../data/' dataFolder dataNameSL '.trc']);
%                     markerNamesSL = fieldnames(trcSL.data);
%                     markerNamesSL = markerNamesSL(3:end); % first 2 names are Frame # and Time
%                     % Rotate markers so subject jumps in positive X direction
%                     for m = 1:length(markerNamesSL)
%                         trcSL.data.(markerNamesSL{m}) = (rotz(-pi/2)*trcSL.data.(markerNamesSL{m})')';
%                     end
%                     trcSL = fillInTRCFrames(trcSL);
% 
%                     %Determine TO and landing frames in recording
                    dt = 1/(trc.DataRate);                
%                     footSpeedX = diff(trc.data.FRMed(:,1));
%                     [~,midJumpFrame] = max(footSpeedX);
%                     footTOFrame = find(footSpeedX(1:midJumpFrame)<=1,1,'last'); % equivalent to (1/1000)/dt = 0.2 [m/s]
%                     footLandFrame = midJumpFrame + find(footSpeedX(midJumpFrame:end)<=1,1,'first');
%                     startFrame = footTOFrame - trc.DataRate; %add 1 second to start
%                     endFrame = footLandFrame + trc.DataRate; %add 1 second to end



                    %Create Euler and Lie Group Models
                    initPosWindow = [101, 200];
%                     [mdl_eul, trcMod_eul, initPos] = createJumpModel(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);
                    [mdl_eul, trcMod_eul, initPos, ~] = createJumpModel_v2(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);

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
                    
                    ekf_eul = EKF_Q_DQ_DDQ(mdl_eul);
%                     ekf_eul = EKF_Q_DQ(mdl_eul);

                    ekf_eul.observation_noise = eye(ekf_eul.sizeZ) * 0.01;

                    ekf_eul.covariance = eye(ekf_eul.sizeX) * 1;

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
                    ekf_eul.process_noise = eye(ekf_eul.sizeX) * x;

                    ekf_covariance = zeros(size(ekf_eul.covariance,1),size(ekf_eul.covariance,2),size(mes_eul,1));

                    %We will also save hatInv of state instead of full out state
                    ekf_state = zeros(size(mes_eul,1),ekf_eul.sizeX);

                    jointAngles_names = {mdl_eul.joints.name};
                    jointAngles_eul = [];


                    % For new marker swapping MatlabWrapper
                    mes_obj = SensorMeasurement(1,30);
                    [mes_obj.size] = deal(3);
                    [mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);
                    matches = 1:numel(mdl_eul.sensors);
                    matches = [matches' matches'];


                    for i=1:size(mes_eul,1)

                        z_eul = mes_eul(i,:);
                        mes_obj.setMesArray(z_eul); % For marker swapping MatlabWrapper
                        u = dt;

                        ekf_eul.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper
                        mes_eul_predict(i,:) = ekf_eul.makeMeasure(ekf_eul.state).getMesArray(); % For marker swapping MatlabWrapper

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


                    currJump = 12*(i_targ-1) + 6*(i_set-1) + i_jump;
                    markerMAE_all_ave(i_part,currJump) = mean(mae_eul_FB);
                    markerMAE_marker_ave(i_part,currJump,:) = mae_eul_markers;


                    disp([dataName ' filtering complete'])
                end
    % END MAIN NESTED FOR LOOPS
            end
        end
    end

    % finished this participant
end
% EKF_allPart_runtime = toc

markerMAE_all_ave(markerMAE_all_ave==0) = NaN;

f = nanmean(nanmean(markerMAE_all_ave));
