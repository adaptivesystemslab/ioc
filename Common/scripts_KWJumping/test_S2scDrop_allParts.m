%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
% clear vis; close all;

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);

visualize = 0;
plotJointAngles = 0;
pauseActive = 0;

S2sc_L_mean = zeros(22,1);
S2sc_R_mean = zeros(22,1);
% S2sc_L3_mean = zeros(22,3);
% S2sc_R3_mean = zeros(22,3);
S2sc_L_var = zeros(22,1);
S2sc_R_var = zeros(22,1);
% S2sc_L3_var = zeros(22,3);
% S2sc_R3_var = zeros(22,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
%               '13','14','15','16','17','18','19','20','21','22'};
partsToRun = {'05','07','09','21'};

% Large Nested FOR Loops for each jump recording
for i_part = 1:numel(partsToRun)
    partNum = partsToRun{i_part};
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};

    data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));
    
    
    S2sc_L = zeros(36,1);
    S2sc_R = zeros(36,1);
%     S2sc_L3 = zeros(36,3);
%     S2sc_R3 = zeros(36,3);
    
    
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
                    
                    currJump = 12*(i_targ-1) + 6*(i_set-1) + i_jump;
                    
                    
                    
%                     window = [50,600];
                    window = [251,450];
                    [S2sc_L(currJump), S2sc_R(currJump), normUpArm_L, normUpArm_R] = getS2scDrop(trc,0,window);
%                     [S2sc_L3(currJump,:), S2sc_R3(currJump,:)] = getS2scDrop(trc,1);
                    
                    figure(1); clf; hold on; grid on;
                    plot(normUpArm_L,'r');
                    plot(normUpArm_R,'b');
                    ylim([0.15, 0.45]);
                    pause(0.01);
                end
            end
        end
        
    end
    
    S2sc_L_mean(i_part+1) = mean(S2sc_L);
    S2sc_R_mean(i_part+1) = mean(S2sc_R);
    S2sc_L_var(i_part+1) = var(S2sc_L);
    S2sc_R_var(i_part+1) = var(S2sc_R);
    
    
%     figure(2); clf;
%     subplot(2,1,1); hold on; grid on;    
%     plot(S2sc_L,'r');
%     plot([1,36],[S2sc_L_mean(i_part+1),S2sc_L_mean(i_part+1)],'k--');
%     ylim([-0.09,-0.02]);
%     title('Left Shoulder');
%     subplot(2,1,2); hold on; grid on;     
%     plot(S2sc_R,'r');
%     ylim([-0.09,-0.02]);
%     title('Right Shoulder');
%     plot([1,36],[S2sc_R_mean(i_part+1),S2sc_R_mean(i_part+1)],'k--');
    
%     w = waitforbuttonpress();
%     pause(0.01);
    
    disp(['P' partNum ' done']);
end




