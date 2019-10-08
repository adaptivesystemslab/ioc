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

partNum = '06';

targetDist = {'55','70','85'};
markerNames = {'Neck';'SL';'SR';'HL';'HR';'BT';'BL';'BR';...
    'ELLat';'ELMed';'WLLat';'WLMed';'ERLat';'ERMed';'WRLat';'WRMed';...
    'KLLat';'KLMed';'ALLat';'ALMed';'FLLat';'FLMed';'FLHeel';...
    'KRLat';'KRMed';'ARLat';'ARMed';'FRLat';'FRMed';'FRHeel'};

mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end));

marker_mae_FB = zeros(36,1);

for i_targ = 3%1:3
    for i_set = 2%1:2
        for i_jump = 3%1:6
            loadFilePath = [newdir 'results\RESULTS_EKF_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
            if(exist(loadFilePath,'file')~=2)
                disp(['RESULTS_EKF_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat NOT FOUND']);
            else
                load(loadFilePath);
                
                
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
                
                
                
                
                %% Compute the Errors
                err_eul = mes_eul - mes_eul_predict;
%                 err_lie = mes_lie - mes_lie_predict;
                
                %Distance between each marker real and prediction, at each timestep
                norm_err_eul = zeros(size(err_eul,1),size(err_eul,2)/3);
%                 norm_err_lie = zeros(size(err_lie,1),size(err_lie,2)/3);
                
                
                for i=1:(size(err_eul,2)/3)
                    norm_err_eul(:,i) = sqrt(sum(err_eul(:,i*3-2:i*3).^2,2));
%                     norm_err_lie(:,i) = sqrt(sum(err_lie(:,i*3-2:i*3).^2,2));
                end
                
                
                % MAE for full body (FB), at each timestep
                mae_eul_FB = zeros(size(norm_err_eul,1),1);
%                 mae_lie_FB = zeros(size(norm_err_lie,1),1);
                for i = 1:size(norm_err_eul,1)
                    mae_eul_FB(i) = mean(norm_err_eul(i,:));
%                     mae_lie_FB(i) = mean(norm_err_lie(i,:));
                end
                
                % MAE for each marker, averaged over all time
                mae_eul_markers = zeros(size(norm_err_eul,2),1);
%                 mae_lie_markers = zeros(size(norm_err_lie,2),1);
                for i = 1:size(norm_err_eul,2)
                    mae_eul_markers(i) = mean(norm_err_eul(:,i));
%                     mae_lie_markers(i) = mean(norm_err_lie(:,i));
                end
                
                
                %% Make Figures
                
                if(plotMAE_allMarkers)
                    %Marker MAE, all 30 markers over time
                    figure(1); clf; hold on;
                    % First EUL Model
%                     subplot(2,1,1); hold on;
                    % torso and neck markers
                    plot(norm_err_eul(:,1),'k'); %neck
                    plot(norm_err_eul(:,4),'k--'); %hips
                    plot(norm_err_eul(:,5),'k--');
                    plot(norm_err_eul(:,6),'k-.'); %back
                    plot(norm_err_eul(:,7),'k-.');
                    plot(norm_err_eul(:,8),'k-.');
                    % left arm markers
                    plot(norm_err_eul(:,2),'b'); %shoulder
                    plot(norm_err_eul(:,9),'b--'); %elbows
                    plot(norm_err_eul(:,10),'b--');
                    plot(norm_err_eul(:,11),'b-.'); %wrists
                    plot(norm_err_eul(:,12),'b-.');
                    % right arm markers
                    plot(norm_err_eul(:,3),'r');
                    plot(norm_err_eul(:,13),'r--');
                    plot(norm_err_eul(:,14),'r--');
                    plot(norm_err_eul(:,15),'r-.');
                    plot(norm_err_eul(:,16),'r-.');
                    % left leg markers
                    plot(norm_err_eul(:,17),'c'); %knees
                    plot(norm_err_eul(:,18),'c');
                    plot(norm_err_eul(:,19),'c--'); %ankles
                    plot(norm_err_eul(:,20),'c--');
                    plot(norm_err_eul(:,21),'c-.'); %feet
                    plot(norm_err_eul(:,22),'c-.');
                    plot(norm_err_eul(:,23),'c-.');
                    % right leg markers
                    plot(norm_err_eul(:,24),'g');
                    plot(norm_err_eul(:,25),'g');
                    plot(norm_err_eul(:,26),'g--');
                    plot(norm_err_eul(:,27),'g--');
                    plot(norm_err_eul(:,28),'g-.');
                    plot(norm_err_eul(:,29),'g-.');
                    plot(norm_err_eul(:,30),'g-.');
                    % TO and Landing frames
                    plot([footTOFrame,footTOFrame],[0,0.5],'k--');
                    plot([footLandFrame,footLandFrame],[0,0.5],'k--');
                    plot([startFrame,startFrame],[0,0.5],'k--');
                    plot([endFrame,endFrame],[0,0.5],'k--');

                    ylabel('Distance [m]');
                    axis([0 size(norm_err_eul,1) 0 0.06]);
                    grid on;


%                     % Second LIE Model
%                     subplot(2,1,2); hold on;
%                     % torso and neck markers
%                     plot(norm_err_lie(:,1),'k');
%                     plot(norm_err_lie(:,4),'k--');
%                     plot(norm_err_lie(:,5),'k--');
%                     plot(norm_err_lie(:,6),'k-.');
%                     plot(norm_err_lie(:,7),'k-.');
%                     plot(norm_err_lie(:,8),'k-.');
%                     % left arm markers
%                     plot(norm_err_lie(:,2),'b');
%                     plot(norm_err_lie(:,9),'b--');
%                     plot(norm_err_lie(:,10),'b--');
%                     plot(norm_err_lie(:,11),'b-.');
%                     plot(norm_err_lie(:,12),'b-.');
%                     % right arm markers
%                     plot(norm_err_lie(:,3),'r');
%                     plot(norm_err_lie(:,13),'r--');
%                     plot(norm_err_lie(:,14),'r--');
%                     plot(norm_err_lie(:,15),'r-.');
%                     plot(norm_err_lie(:,16),'r-.');
%                     % left leg markers
%                     plot(norm_err_lie(:,17),'c');
%                     plot(norm_err_lie(:,18),'c');
%                     plot(norm_err_lie(:,19),'c--');
%                     plot(norm_err_lie(:,20),'c--');
%                     plot(norm_err_lie(:,21),'c-.');
%                     plot(norm_err_lie(:,22),'c-.');
%                     plot(norm_err_lie(:,23),'c-.');
%                     % right leg markers
%                     plot(norm_err_lie(:,24),'g');
%                     plot(norm_err_lie(:,25),'g');
%                     plot(norm_err_lie(:,26),'g--');
%                     plot(norm_err_lie(:,27),'g--');
%                     plot(norm_err_lie(:,28),'g-.');
%                     plot(norm_err_lie(:,29),'g-.');
%                     plot(norm_err_lie(:,30),'g-.');
%                     % TO and Landing frames
%                     plot([footTOFrame,footTOFrame],[0,0.5],'k--');
%                     plot([footLandFrame,footLandFrame],[0,0.5],'k--');
%                     plot([startFrame,startFrame],[0,0.5],'k--');
%                     plot([endFrame,endFrame],[0,0.5],'k--');

                    suptitle(['P' partNum '-target' num2str(i_targ) '-jump' num2str(6*(i_set-1)+i_jump) ': MAE Of All Markers']);
                    xlabel('Timestep');
                    ylabel('Distance [m]');
                    xlim([0 size(norm_err_eul,1)]);
                    grid on;
                    
                    if(pauseActive)
                        w = waitforbuttonpress();
                    end
                end
                
                
%                 markerDiff = norm_err_eul(startFrame:endFrame,:)-norm_err_lie(startFrame:endFrame,:);
%                 [maxDiff_frame,maxDiff_marker] = find(max(max(abs(markerDiff)))==abs(markerDiff));
%                 maxDiff_frame = maxDiff_frame + startFrame;
%                 disp(['P' partNum '_' num2str(i_targ) '_' num2str(6*(i_set-1)+i_jump)...
%                             ' Max Marker Diff [mm]: ' num2str(max(max(abs(markerDiff)))*1000)...
%                             ', by marker ' markerNames{maxDiff_marker} ', at frame ' num2str(maxDiff_frame)])
                
                
                if(plotMAE_markersPerLimb)
                    
                    limbList = {'torso/neck','larm','rarm','lleg','rleg'};
                    markerList = {[1,4:8],[2,9:12],[3,13:16],[17:23],[24:30]};
                    
                    figure(2); 
                    for i = 1:5
                        clf; hold on;
                        for j = 1:numel(markerList{i})
                            plot(norm_err_eul(:,markerList{i}(j)),'r');
                        end
%                         for j = 1:numel(markerList{i}) %eul and lie plot commands separated to better format legend
%                             plot(norm_err_lie(:,markerList{i}(j)),'b');
%                         end
                        
                        plot([footTOFrame,footTOFrame],[0,0.1],'k--');
                        plot([footLandFrame,footLandFrame],[0,0.1],'k--');
                        plot([startFrame,startFrame],[0,0.1],'k--');
                        plot([endFrame,endFrame],[0,0.1],'k--');
                        
                        title(['P' partNum '-target' num2str(i_targ) '-jump' num2str(6*(i_set-1)+i_jump) ': ' limbList{i}]);
                        legend({markerNames{markerList{i}},'blue=Lie'});
                        xlabel('Timestep');
                        ylabel('Distance [m]');
                        axis([0 size(norm_err_eul,1) 0 0.05]);
                        grid on;
                        
                        
                        maxMarkerDiff_pos = max(max(markerDiff(:,markerList{i})));
                        maxMarkerDiff_neg = min(min(markerDiff(:,markerList{i})));
                        disp(['P' partNum '-target' num2str(i_targ) '-jump' num2str(6*(i_set-1)+i_jump)...
                            ' Max Marker Diff [mm]: ' num2str(maxMarkerDiff_pos*1000) ', ' num2str(maxMarkerDiff_neg*1000)])
                
                        
                        if(pauseActive)
                            w=0;
                            while(~w)
                                w = waitforbuttonpress(); %waits for keyboard press
                            end
                        end
                    end
                    
                    
                end
                
                
                
                if(plotMAE_FB)
                    %Marker MAE, full body
                    figure(3); clf; hold on; grid on;
                    plot(mae_eul_FB);
%                     plot(mae_lie_FB);
                    % TO and Landing frames
                    plot([footTOFrame,footTOFrame],[0,0.05],'k--');
                    plot([footLandFrame,footLandFrame],[0,0.05],'k--');
                    plot([startFrame,startFrame],[0,0.05],'k--');
                    plot([endFrame,endFrame],[0,0.05],'k--');
                    
%                     legend('EKF','LG EKF');
                    title(['P' partNum '-target' num2str(i_targ) '-jump' num2str(6*(i_set-1)+i_jump) ': MAE Of All Markers']);
                    xlabel('Timestep');
                    ylabel('Distance [m]');
                    axis([0, size(mae_eul_FB,1), 0, 0.05]);
                    grid on;
                    
                    if(pauseActive)
                        w = waitforbuttonpress();
                    end
                end
                
                
%                 %Marker MAE, torso and each limb
%                 figure(3); clf;
%                 titles = {'rarm','larm','torso/head','rleg','lleg'};
%                 for i = 1:5
%                     subplot(3,2,i); hold on; grid on;
%                     plot(mae_eul_limbs(:,i));
% %                     hold on
% %                     plot(mae_lie_limbs(:,i));
% %                     legend('EKF','LG EKF');
%                     title([titles{i} 'Markers MAE' num2str(1/dt) ' Hz)']);
%                     ylabel('Distance [m]');
%                     xlim([0 size(mae_eul_limbs,1)]);
%                 end
%                 xlabel('Timestep');
                
                
                
                % Mean MAE
                
                disp('Mean MAE values [mm] = ')
                disp(['EKF, ', num2str(1000*mean(mae_eul_FB))])
%                 disp(['LG EKF, ', num2str(1000*mean(mae_lie_FB(5:end)))])
%                 mae_diff = 1000*(mean(mae_eul_FB(5:end)) - mean(mae_lie_FB(5:end)));
%                 disp(['Eul-LG mae diff [mm], ', num2str(mae_diff) ])
                jumpNum = (i_targ-1)*12 + (i_set-1)*6 + i_jump;
                marker_mae_FB(jumpNum,1) = mean(mae_eul_FB(startFrame:endFrame,:));
                
                pause(0.1);
            end
            
            
            
            
        end
    end
end


marker_mae_FB;



