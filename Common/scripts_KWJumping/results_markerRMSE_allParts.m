% Reads in previously-saved trc and ekf marker position data, calculates
% mae error and makes plots

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
addpath('..\..\');
addpath(EKFCodePath);

% Choose dataset and parameters
plotMAE_allMarkers = 1;
plotMAE_markersPerLimb = 0;
plotMAE_FB = 0;
pauseActive = 1;

partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'16'};



mae_ave = zeros(22,30);
mae_ave_jumpOnly = zeros(22,30);
mae_ave_constLinks = zeros(22,30);
mae_ave_jumpOnly_constLinks = zeros(22,30);

mae_max_constLinks = zeros(22,30);


for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};

%     markerNames = {'Neck';'SL';'SR';'HL';'HR';'BT';'BL';'BR';...
%         'ELLat';'ELMed';'WLLat';'WLMed';'ERLat';'ERMed';'WRLat';'WRMed';...
%         'KLLat';'KLMed';'ALLat';'ALMed';'FLLat';'FLMed';'FLHeel';...
%         'KRLat';'KRMed';'ARLat';'ARMed';'FRLat';'FRMed';'FRHeel'};

    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end));
    
    
    markerMAEPerJump = zeros(36,30);
    markerMAEPerJump_constLinks = zeros(36,30);
    
    for i_targ = 1:3
        for i_set = 1:2
            for i_jump = 1:6
                loadFilePath = [newdir 'results\RESULTS_EKF_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
                loadFilePath2 = [newdir 'results\RESULTS_EKF_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '_const_model.mat'];
                if(exist(loadFilePath,'file')~=2)
                    disp(['RESULTS_EKF_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat NOT FOUND']);
                else
                    %% Read in marker data to calculate TO and Landing frames
                    dataFolder = ['P' partNum '_filtered/'];
                    dataName = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                            '_' num2str(i_jump) '_clean-P' partNum '_template'];
                    dataNameSL = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                            '_' num2str(i_jump) '_clean-start_land'];

                    trc = readTrc(['../data/' dataFolder dataName '.trc']);
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
    %                 trcSL = fillInTRCFrames(trcSL);

                    %Determine TO and landing frames in recording
                    dt = 1/(trc.DataRate);
                    footSpeedX = diff(trc.data.FRMed(:,1));
                    [~,midJumpFrame] = max(footSpeedX);
                    footTOFrame = find(footSpeedX(1:midJumpFrame)<=1,1,'last'); % equivalent to (1/1000)/dt = 0.2 [m/s]
                    footLandFrame = midJumpFrame + find(footSpeedX(midJumpFrame:end)<=1,1,'first');
                    startFrame = footTOFrame - trc.DataRate; %add 1 second to start
                    endFrame = footLandFrame + trc.DataRate; %add 1 second to end
                    if(endFrame > size(trc.data.FRMed,1))
                        endFrame = size(trc.data.FRMed,1);
                    end

                    
                    %% load EKF data for changing model
                    load(loadFilePath);

                    %% Compute the Errors
                    err_eul = mes_eul - mes_eul_predict;
    %                 err_lie = mes_lie - mes_lie_predict;

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
                    mae_eul_markers = zeros(size(norm_err_eul,2)-5,1);
                    for i = 1:size(norm_err_eul,2)
                        mae_eul_markers(i) = mean(norm_err_eul(6:end,i));
                    end



                    % Mean MAE

%                     disp('Mean MAE values [mm] = ')
%                     disp(['EKF, ', num2str(1000*mean(mae_eul_FB))])
    %                 disp(['LG EKF, ', num2str(1000*mean(mae_lie_FB(5:end)))])
    %                 mae_diff = 1000*(mean(mae_eul_FB(5:end)) - mean(mae_lie_FB(5:end)));
    %                 disp(['Eul-LG mae diff [mm], ', num2str(mae_diff) ])
                    currJump = (i_targ-1)*12 + (i_set-1)*6 + i_jump;
%                     mae_ave(str2num(partNum),currJump) = mean(mae_eul_FB);
%                     mae_ave_jumpOnly(str2num(partNum),currJump) = mean(mae_eul_FB(startFrame:endFrame,:));
                    markerMAEPerJump(currJump,:) = mae_eul_markers';
                    
                    
                    %% load EKF data for constant link length model
                    load(loadFilePath2);
                    
                     %% Compute the Errors
                    err_eul = mes_eul - mes_eul_predict;
    %                 err_lie = mes_lie - mes_lie_predict;

                    %Distance between each marker real and prediction, at each timestep
                    norm_err_eul = zeros(size(err_eul,1),size(err_eul,2)/3);

                    for i=1:(size(err_eul,2)/3)
                        norm_err_eul(:,i) = sqrt(sum(err_eul(:,i*3-2:i*3).^2,2));
                    end
                    plot(norm_err_eul(:,2:3)); % shoulder error

                    % MAE for full body (FB), at each timestep
                    mae_eul_FB = zeros(size(norm_err_eul,1),1);
                    for i = 1:size(norm_err_eul,1)
                        mae_eul_FB(i) = mean(norm_err_eul(i,:));
                    end

                    % MAE for each marker, averaged over all time
                    mae_eul_markers = zeros(size(norm_err_eul,2)-5,1);
                    for i = 1:size(norm_err_eul,2)
                        mae_eul_markers(i) = mean(norm_err_eul(6:end,i));
                    end
                    
                    
%                     mae_ave_constLinks(str2num(partNum),currJump) = mean(mae_eul_FB);
%                     mae_ave_jumpOnly_constLinks(str2num(partNum),currJump) = mean(mae_eul_FB(startFrame:endFrame,:));
                    markerMAEPerJump_constLinks(currJump,:) = mae_eul_markers';
                    
                    
                    % Max MAE for each marker, throughout single jump
                    mae_max_markers = zeros(size(norm_err_eul,2)-5,1);
                    for i = 1:size(norm_err_eul,2)
                        mae_max_markers(i) = max(norm_err_eul(startFrame:endFrame,i));
                    end
                    
                    marker_maxMAE_constLinks(currJump,:) = mae_max_markers';
                    
%                     figure(1);
%                     for m = 1:30
%                         clf; hold on; grid on;
%                         plot(norm_err_eul(:,m));
%                         title([partNum ' ' markerNames{m}]);
%                         ylim([0 0.06]);
%                         w = waitforbuttonpress();
%                     end
                    
%                     figure(2); clf; hold on; grid on;
%                     plot([repmat(mae_eul_FB(5),5,1); mae_eul_FB(5:end)],'r','LineWidth',2);
%                     plot([footTOFrame,footTOFrame],[0,0.02],'g--','LineWidth',2);
%                     plot([footLandFrame,footLandFrame],[0,0.02],'k--','LineWidth',2);
%                     title(['P' partNum ' Full Body Average MAE']);
%                     xlabel('Frame');
%                     ylabel('Marker MAE [m]');
%                     ylim([0 0.02]);
                    
                    
                end




            end
        end
    end
    
    
    mae_ave(str2num(partNum),:) = mean(markerMAEPerJump,1);
    mae_ave_constLinks(str2num(partNum),:) = mean(markerMAEPerJump_constLinks,1);
    
    mae_max_constLinks(str2num(partNum),:) = mean(marker_maxMAE_constLinks);
    
    disp(['Part ' partNum ' done']);
    
end


% Any parts not calculated still have rows of zeros
mae_ave(mae_ave==0) = NaN;
mae_ave_constLinks(mae_ave_constLinks==0) = NaN;
% mae_ave_jumpOnly(mae_ave_jumpOnly==0) = NaN;
% mae_ave_jumpOnly_constLinks(mae_ave_jumpOnly_constLinks==0) = NaN;
mae_max_constLinks(mae_max_constLinks==0) = NaN;


% Average results for all frames
meanImproveOfConstModel = nanmean(nanmean(mae_ave - mae_ave_constLinks))
maxImproveOfConstModel = max(max(mae_ave - mae_ave_constLinks))
minImproveOfConstModel = min(min(mae_ave - mae_ave_constLinks))
mae_diff = mae_ave - mae_ave_constLinks;

% Average results for jump frames only (between startFrame and endFrame)
% meanImproveOfConstModel_jumpOnly = nanmean(nanmean(mae_ave_jumpOnly - mae_ave_jumpOnly_constLinks))
% maxImproveOfConstModel_jumpOnly = max(max(mae_ave_jumpOnly - mae_ave_jumpOnly_constLinks))
% minImproveOfConstModel_jumpOnly = min(min(mae_ave_jumpOnly - mae_ave_jumpOnly_constLinks))
