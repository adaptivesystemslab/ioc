% Reads in previously-saved joint angle data and makes plots, compares
% joint angles of one participant over consecutive jumps (see "learning")
clear;

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

% Choose dataset and parameters

plotJointAngles = 1;
% fixOutliers = 0; % SHOULDN'T NEED THIS IF MARKER DATA IS PROPERLY LABELLED AND CLEANED
saveJAStruct = 0;
use_const_model_links = 1;
use_shldrPrismModel = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
%               '13','14','15','16','17','18','19','20','21','22'};
partsToRun = {'04'};

% flightMean = zeros(numel(partsToRun),4);
% flightVar = zeros(numel(partsToRun),4);

mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 


for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
    targetDist = {'55','70','85'};

    JA = []; %"joint angles all"
    JA.partNum = partNum;
    JA.missingData = zeros(12,3);
%     JA.noisyData = zeros(12,3);
%     JA.noisyData_lie = zeros(12,3);
%     JA.shiftFlipData = zeros(36,33);
%     JA.shiftFlipData_lie = zeros(36,33);
    numJoints = 0;
    modelLinkNames = [];

    %% Read in and set participant data and jump grading
    [age,gender,height,weight,S2sc] = getPartData(partNum);
    JA.partData.age = age;
    JA.partData.gender = gender;
    JA.partData.height = double(height)/100; %to [m]
    JA.partData.weight = weight;
    JA.partData.S2sc = S2sc/100; %to [m]
    [targetLengths,jumpGrades] = getJumpGrading(partNum);
    JA.targetLengths = targetLengths;
    JA.jumpGrades = jumpGrades;
    
    JA.modelLinks = [];
    JA.world2base = [];
    
    
    %% Read in joint angle data for all jumps for this participant
    for i_targ = 1:3
        for i_set = 1:2
            for i_jump = 1:6
                currJump = 6*(i_set-1) + i_jump;

                
                saveFilePath = [newdir '\results\RESULTS_JA_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
                if(exist(saveFilePath,'file')==2)
                    load(saveFilePath);

                    JA.targ(i_targ).jump(currJump).data = jointAngles_eul;
    %                 JA.targ_lie(i_targ).jump(currJump).data = jointAngles_lie;
                    JA.dataLengths(currJump,i_targ) = size(jointAngles_eul,1);
                    JA.startFrame(currJump,i_targ) = startFrame;
                    JA.TOFrame(currJump,i_targ) = footTOFrame;
                    JA.LandFrame(currJump,i_targ) = footLandFrame;
                    JA.endFrame(currJump,i_targ) = endFrame;
                    
                    % rotation matricies for hips and shoulders
                    JA.targ_R_mat(i_targ).jump(currJump) = R_mat;
                    
                    if(isempty(modelLinkNames))
                        if(use_const_model_links) % only need to do this once
                            load([newdir '\results\Model_Info_JA_P' partNum '.mat']);
                            modelLinkNames = fieldnames(modelLinks_ave);
                            
                            for link = 1:numel(modelLinkNames)
                                JA.modelLinks.(modelLinkNames{link}) = modelLinks_ave.(modelLinkNames{link});
                            end
                        end
                    end
                    
                    if(~use_const_model_links)
                        modelLinkNames = fieldnames(modelLinks);
                        
                        for link = 1:numel(modelLinkNames)
                            JA.modelLinks.(modelLinkNames{link})(12*(i_targ-1)+currJump,:) = ...
                                modelLinks.(modelLinkNames{link});
                        end
                    end
                    
                    JA.world2base(12*(i_targ-1)+currJump,:,:) = world2base;
                    JA.locationStart(12*(i_targ-1)+currJump) = mean(locationStart(:,1))'; % global x-value of starting line
                    JA.locationLand(12*(i_targ-1)+currJump,:) = mean(locationLand(:,1:2));
                    
                    if(numJoints==0)
                        JA.jointNames = jointAngles_names;
                        numJoints = size(jointAngles_eul,2);
                    end

                else %file not available, populate JAall with zeros
                    JA.missingData(6*(i_set-1) + i_jump,i_targ) = 1;

                    if(numJoints==0)    numJoints = 33;     end %first file missing, assume 33 joints
                    JA.targ(i_targ).jump(currJump).data = zeros(2500,numJoints);
    %                 JA.targ_lie(i_targ).jump(currJump).data = zeros(2500,numJoints);
                    JA.dataLengths(currJump,i_targ) = 2500;
                    JA.startFrame(currJump,i_targ) = 1300;
                    JA.TOFrame(currJump,i_targ) = 1500;
                    JA.LandFrame(currJump,i_targ) = 1600;
                    JA.endFrame(currJump,i_targ) = 1800;
                end
            end
        end
    end
    
    if(strcmp(partNum,'09'))
        JA.locationStart(25) = JA.locationStart(29); % start line markers not moved from T3 S2 until 5th jump of set
        JA.locationStart(26) = JA.locationStart(29);
        JA.locationStart(27) = JA.locationStart(29);
        JA.locationStart(28) = JA.locationStart(29);
    end
    
%     flightTime = JA.LandFrame - JA.TOFrame;
%     flightMean(i_part,1:3) = mean(flightTime);
%     flightVar(i_part,1:3) = var(flightTime);
%     disp(['Part. ' partNum]);
%     if(max(max(flightTime)) > 130)
%         flightTime
%         disp(['Part. ' partNum]);
%     end
    
    
    %% Participant 2 uses T-pose calibration, and gimbal lock messes up shoulders, manually fix here
    if(partNum == '02')
        if(~use_shldrPrismModel)
            % Target 1
            % rshldrElev
            JA.targ(1).jump(1).data(:,10) = JA.targ(1).jump(1).data(:,10) + pi;
            % lshldrElev
            JA.targ(1).jump(1).data(:,15) = JA.targ(1).jump(1).data(:,15) + pi;
            JA.targ(1).jump(2).data(:,15) = JA.targ(1).jump(2).data(:,15) + pi;
            JA.targ(1).jump(5).data(:,15) = JA.targ(1).jump(5).data(:,15) + pi;
            JA.targ(1).jump(6).data(:,15) = JA.targ(1).jump(6).data(:,15) + pi;
            JA.targ(1).jump(7).data(:,15) = JA.targ(1).jump(7).data(:,15) + pi;
            JA.targ(1).jump(8).data(:,15) = JA.targ(1).jump(8).data(:,15) + pi;
            % rshldrAbd
            JA.targ(1).jump(1).data(:,11) = -(JA.targ(1).jump(1).data(:,11) - pi);
            % lshldrAbd
            JA.targ(1).jump(1).data(:,16) = -(JA.targ(1).jump(1).data(:,16) - pi);
            JA.targ(1).jump(2).data(:,16) = -(JA.targ(1).jump(2).data(:,16) - pi);
            JA.targ(1).jump(5).data(:,16) = -(JA.targ(1).jump(5).data(:,16) - pi);
            JA.targ(1).jump(6).data(:,16) = -(JA.targ(1).jump(6).data(:,16) - pi);
            JA.targ(1).jump(7).data(:,16) = -(JA.targ(1).jump(7).data(:,16) - pi);
            JA.targ(1).jump(8).data(:,16) = -(JA.targ(1).jump(8).data(:,16) - pi);
            % rshldrExtRot
            JA.targ(1).jump(1).data(:,12) = JA.targ(1).jump(1).data(:,12) - pi;
            % lshldrExtRot
            JA.targ(1).jump(1).data(:,17) = JA.targ(1).jump(1).data(:,17) - pi;
            JA.targ(1).jump(2).data(:,17) = JA.targ(1).jump(2).data(:,17) - pi;
            JA.targ(1).jump(5).data(:,17) = JA.targ(1).jump(5).data(:,17) - pi;
            JA.targ(1).jump(6).data(:,17) = JA.targ(1).jump(6).data(:,17) - pi;
            JA.targ(1).jump(7).data(:,17) = JA.targ(1).jump(7).data(:,17) - pi;
            JA.targ(1).jump(8).data(:,17) = JA.targ(1).jump(8).data(:,17) - pi;

            % Target 2
            % lshldrElev
            JA.targ(2).jump(1).data(:,15) = JA.targ(2).jump(1).data(:,15) + pi;
            JA.targ(2).jump(3).data(:,15) = JA.targ(2).jump(3).data(:,15) + pi;
            JA.targ(2).jump(12).data(:,15) = JA.targ(2).jump(12).data(:,15) + pi;
            % lshldrAbd
            JA.targ(2).jump(1).data(:,16) = -(JA.targ(2).jump(1).data(:,16) - pi);
            JA.targ(2).jump(3).data(:,16) = -(JA.targ(2).jump(3).data(:,16) - pi);
            JA.targ(2).jump(12).data(:,16) = -(JA.targ(2).jump(12).data(:,16) - pi);
            % lshldrExtRot
            JA.targ(2).jump(1).data(:,17) = JA.targ(2).jump(1).data(:,17) - pi;
            JA.targ(2).jump(3).data(:,17) = JA.targ(2).jump(3).data(:,17) - pi;
            JA.targ(2).jump(12).data(:,17) = JA.targ(2).jump(12).data(:,17) - pi;

            % Target 3
            % lshldrElev
            JA.targ(3).jump(1).data(:,15) = JA.targ(3).jump(1).data(:,15) + pi;
            JA.targ(3).jump(2).data(:,15) = JA.targ(3).jump(2).data(:,15) + pi;
            JA.targ(3).jump(6).data(:,15) = JA.targ(3).jump(6).data(:,15) - pi;
            JA.targ(3).jump(7).data(:,15) = JA.targ(3).jump(7).data(:,15) + pi;
            % lshldrAbd
            JA.targ(3).jump(1).data(:,16) = -(JA.targ(3).jump(1).data(:,16) - pi);
            JA.targ(3).jump(2).data(:,16) = -(JA.targ(3).jump(2).data(:,16) - pi);
            JA.targ(3).jump(6).data(:,16) = -(JA.targ(3).jump(6).data(:,16) - pi);
            JA.targ(3).jump(7).data(:,16) = -(JA.targ(3).jump(7).data(:,16) - pi);
            % lshldrExtRot
            JA.targ(3).jump(1).data(:,17) = JA.targ(3).jump(1).data(:,17) - pi;
            JA.targ(3).jump(2).data(:,17) = JA.targ(3).jump(2).data(:,17) - pi;
            JA.targ(3).jump(6).data(:,17) = JA.targ(3).jump(6).data(:,17) + pi;
            JA.targ(3).jump(7).data(:,17) = JA.targ(3).jump(7).data(:,17) - pi;
        else %using shldrPrism model
            % Target 1
            % rshldrElev
            JA.targ(1).jump(1).data(:,11) = JA.targ(1).jump(1).data(:,10) + pi;
            % lshldrElev
            JA.targ(1).jump(1).data(:,17) = JA.targ(1).jump(1).data(:,15) + pi;
            JA.targ(1).jump(2).data(:,17) = JA.targ(1).jump(2).data(:,15) + pi;
            JA.targ(1).jump(5).data(:,17) = JA.targ(1).jump(5).data(:,15) + pi;
            JA.targ(1).jump(6).data(:,17) = JA.targ(1).jump(6).data(:,15) + pi;
            JA.targ(1).jump(7).data(:,17) = JA.targ(1).jump(7).data(:,15) + pi;
            JA.targ(1).jump(8).data(:,17) = JA.targ(1).jump(8).data(:,15) + pi;
            % rshldrAbd
            JA.targ(1).jump(1).data(:,12) = -(JA.targ(1).jump(1).data(:,11) - pi);
            % lshldrAbd
            JA.targ(1).jump(1).data(:,18) = -(JA.targ(1).jump(1).data(:,16) - pi);
            JA.targ(1).jump(2).data(:,18) = -(JA.targ(1).jump(2).data(:,16) - pi);
            JA.targ(1).jump(5).data(:,18) = -(JA.targ(1).jump(5).data(:,16) - pi);
            JA.targ(1).jump(6).data(:,18) = -(JA.targ(1).jump(6).data(:,16) - pi);
            JA.targ(1).jump(7).data(:,18) = -(JA.targ(1).jump(7).data(:,16) - pi);
            JA.targ(1).jump(8).data(:,18) = -(JA.targ(1).jump(8).data(:,16) - pi);
            % rshldrExtRot
            JA.targ(1).jump(1).data(:,13) = JA.targ(1).jump(1).data(:,12) - pi;
            % lshldrExtRot
            JA.targ(1).jump(1).data(:,19) = JA.targ(1).jump(1).data(:,17) - pi;
            JA.targ(1).jump(2).data(:,19) = JA.targ(1).jump(2).data(:,17) - pi;
            JA.targ(1).jump(5).data(:,19) = JA.targ(1).jump(5).data(:,17) - pi;
            JA.targ(1).jump(6).data(:,19) = JA.targ(1).jump(6).data(:,17) - pi;
            JA.targ(1).jump(7).data(:,19) = JA.targ(1).jump(7).data(:,17) - pi;
            JA.targ(1).jump(8).data(:,19) = JA.targ(1).jump(8).data(:,17) - pi;

            % Target 2
            % lshldrElev
            JA.targ(2).jump(1).data(:,17) = JA.targ(2).jump(1).data(:,15) + pi;
            JA.targ(2).jump(3).data(:,17) = JA.targ(2).jump(3).data(:,15) + pi;
            JA.targ(2).jump(12).data(:,17) = JA.targ(2).jump(12).data(:,15) + pi;
            % lshldrAbd
            JA.targ(2).jump(1).data(:,18) = -(JA.targ(2).jump(1).data(:,16) - pi);
            JA.targ(2).jump(3).data(:,18) = -(JA.targ(2).jump(3).data(:,16) - pi);
            JA.targ(2).jump(12).data(:,18) = -(JA.targ(2).jump(12).data(:,16) - pi);
            % lshldrExtRot
            JA.targ(2).jump(1).data(:,19) = JA.targ(2).jump(1).data(:,17) - pi;
            JA.targ(2).jump(3).data(:,19) = JA.targ(2).jump(3).data(:,17) - pi;
            JA.targ(2).jump(12).data(:,19) = JA.targ(2).jump(12).data(:,17) - pi;

            % Target 3
            % lshldrElev
            JA.targ(3).jump(1).data(:,17) = JA.targ(3).jump(1).data(:,15) + pi;
            JA.targ(3).jump(2).data(:,17) = JA.targ(3).jump(2).data(:,15) + pi;
            JA.targ(3).jump(6).data(:,17) = JA.targ(3).jump(6).data(:,15) - pi;
            JA.targ(3).jump(7).data(:,17) = JA.targ(3).jump(7).data(:,15) + pi;
            % lshldrAbd
            JA.targ(3).jump(1).data(:,18) = -(JA.targ(3).jump(1).data(:,16) - pi);
            JA.targ(3).jump(2).data(:,18) = -(JA.targ(3).jump(2).data(:,16) - pi);
            JA.targ(3).jump(6).data(:,18) = -(JA.targ(3).jump(6).data(:,16) - pi);
            JA.targ(3).jump(7).data(:,18) = -(JA.targ(3).jump(7).data(:,16) - pi);
            % lshldrExtRot
            JA.targ(3).jump(1).data(:,19) = JA.targ(3).jump(1).data(:,17) - pi;
            JA.targ(3).jump(2).data(:,19) = JA.targ(3).jump(2).data(:,17) - pi;
            JA.targ(3).jump(6).data(:,19) = JA.targ(3).jump(6).data(:,17) + pi;
            JA.targ(3).jump(7).data(:,19) = JA.targ(3).jump(7).data(:,17) - pi;
        end
        
    end
        



    %% ROUGH ALIGN: Align JAs to land frame, crop frames near start and end of each recording
    JA.targAlign = JA.targ;
    % JA.targAlign_lie = JA.targ_lie;
    dataLengthsAlign = JA.dataLengths;

    for i_targ = 1:numel(JA.targAlign)
        % crop frames off front to align "footLandFrame" for all 12 jumps
        [minLandFrame,idx_LandFrame] = min(JA.LandFrame(:,i_targ));
        JA.LandAlignJump(i_targ) = idx_LandFrame;

        for currJump = 1:numel(JA.targAlign(i_targ).jump)
            if(currJump~=idx_LandFrame)
                cropFrames = JA.LandFrame(currJump,i_targ) - minLandFrame;
                dataLengthsAlign(currJump,i_targ) = dataLengthsAlign(currJump,i_targ) - cropFrames;

                JA.targAlign(i_targ).jump(currJump).data = ...
                    JA.targAlign(i_targ).jump(currJump).data(cropFrames+1:end,:);
    %             JA.targAlign_lie(i_targ).jump(currJump).data = ...
    %                 JA.targAlign_lie(i_targ).jump(currJump).data(cropFrames+1:end,:);
            end
        end

        % crop frames off back to make all 12 jumps the same length
        [minDataLength,idx_dataLength] = min(dataLengthsAlign(:,i_targ));
        for currJump = 1:12
            if(currJump~=idx_dataLength)
                JA.targAlign(i_targ).jump(currJump).data = ...
                    JA.targAlign(i_targ).jump(currJump).data(1:minDataLength,:);
    %             JA.targAlign_lie(i_targ).jump(currJump).data = ...
    %                 JA.targAlign_lie(i_targ).jump(currJump).data(1:minDataLength,:);
            end
        end
    end


    %% FINE ALIGN: temporal shift alignment based on best match of lower body JAs
    [shiftValues, ~] = JAShiftValues(JA); % finds shift values relative to 12th (last) jump to each target
    JA.LandFrame_shiftValues = shiftValues;
    
    for i_targ = 1:numel(JA.targAlign)
        for i = 1:numel(JA.targAlign(i_targ).jump)-1 %use 12th (last) jump as reference
            sig = JA.targAlign(i_targ).jump(i).data;
            shift = shiftValues(i,i_targ);
            if(shift < 0) % shift back in time relative to reference jump
                sig = [zeros(-shift,numJoints); sig(1:end+shift,:)];
            elseif(shift > 0) % shift forward in time relative to reference jump
                sig = [sig(shift+1:end,:); zeros(shift,numJoints)];
            end

            JA.targAlign(i_targ).jump(i).data = sig;
        end
    end


% Below is no longer necessary, fixed all shifting/flipping with mocap
% cleaning and shoulder rotation decomposition. For Participant 02, have
% code block above to manually correct T-pose gimbal lock flipping.

%     %% Vertically shift/flip data, find and remove outliers remove noisy data
%     if(fixOutliers)
%         for i_targ = 1:numel(JA.targAlign)
%             idx_LandFrame = JA.LandAlignJump(i_targ);
%             % get all jump data for specific joint
%             for j = 1:numJoints
%                 sig_eul = zeros(size(JA.targAlign(i_targ).jump(1).data,1),numel(JA.targAlign(i_targ).jump));
%         %         sig_lie = zeros(size(JA.targAlign_lie(i_targ).jump(1).data,1),numel(JA.targAlign_lie(i_targ).jump));
%                 for i = 1:size(sig_eul,2)
%                     sig_eul(:,i) = JA.targAlign(i_targ).jump(i).data(:,j);
%         %             sig_lie(:,i) = JA.targAlign_lie(i_targ).jump(i).data(:,j);
%                 end
% 
%                 % vertically shift/flip JA signals, remove outlier signals, and remove noisy signals
%                 alignWindow = [JA.startFrame(idx_LandFrame,i_targ),JA.endFrame(idx_LandFrame,i_targ)];
%                 [~,shiftSig_eul, shift2Sig_eul, flipSig_eul, shiftFlipSig_eul, shiftFlip2Sig_eul, OL_smooth_eul] ...
%                     = cleanAndAlignData(sig_eul,alignWindow,1);
%     %             [~,shiftSig_lie, shift2Sig_lie, flipSig_lie, shiftFlipSig_lie, shiftFlip2Sig_lie, OL_smooth_lie] ...
%     %                 = cleanAndAlignData(sig_lie,alignWindow,1);
% 
%                 % fix signal shifting/outlier/noise removal for EUL data
%                 for i = 1:numel(JA.targAlign(i_targ).jump)
%                     tmpData = JA.targAlign(i_targ).jump(i).data(:,j);
%                     if(shift2Sig_eul(i)~=0)
%                         JA.targAlign(i_targ).jump(i).data(:,j) = tmpData + 2*pi*shift2Sig_eul(i);
%                         JA.shiftFlipData(12*(i_targ-1) + i,j) = 1;
%                     elseif(shiftSig_eul(i)~=0)
%                         JA.targAlign(i_targ).jump(i).data(:,j) = tmpData + pi*shiftSig_eul(i);
%                         JA.shiftFlipData(12*(i_targ-1) + i,j) = 1;
%                     elseif(shiftFlip2Sig_eul(i)~=0)
%                         JA.targAlign(i_targ).jump(i).data(:,j) = -(tmpData + 2*pi*shiftFlip2Sig_eul(i));
%                         JA.shiftFlipData(12*(i_targ-1) + i,j) = 1;
%                     elseif(shiftFlipSig_eul(i)~=0)
%                         JA.targAlign(i_targ).jump(i).data(:,j) = -(tmpData + pi*shiftFlipSig_eul(i));
%                         JA.shiftFlipData(12*(i_targ-1) + i,j) = 1;
%                     elseif(flipSig_eul(i)~=0)
%                         JA.targAlign(i_targ).jump(i).data(:,j) = -tmpData;
%                         JA.shiftFlipData(12*(i_targ-1) + i,j) = 1;
%                     end
%                     % NOTE: ordered signal mods above from least to most likely to have false positives
% 
%                     if(OL_smooth_eul(i)~=0)
%         %                 JA.targAlign(i_targ).jump(i).data(:,j) = zeros(size(tmpData));
%                         JA.noisyData(i,i_targ) = 1; 
%                     end
%                 end
% 
%         %         % fix signal shifting/outlier/noise removal for LIE data
%         %         for i = 1:numel(JA.targAlign_lie(i_targ).jump)
%         %             tmpData = JA.targAlign_lie(i_targ).jump(i).data(:,j);
%         %             if(shift2Sig_lie(i)~=0)
%         %                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = tmpData + 2*pi*shift2Sig_lie(i);
%         %                 JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
%         %             elseif(shiftSig_lie(i)~=0)
%         %                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = tmpData + pi*shiftSig_lie(i);
%         %                 JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
%         %             elseif(shiftFlip2Sig_lie(i)~=0)
%         %                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = -(tmpData + 2*pi*shiftFlip2Sig_lie(i));
%         %                 JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
%         %             elseif(shiftFlipSig_lie(i)~=0)
%         %                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = -(tmpData + pi*shiftFlipSig_lie(i));
%         %                 JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
%         %             elseif(flipSig_lie(i)~=0)
%         %                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = -tmpData;
%         %                 JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
%         %             end
%         %             % NOTE: ordered signal mods above from least to most likely to have false positives
%         %             
%         %             if(OL_smooth_lie(i)~=0)
%         % %                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = zeros(size(tmpData));
%         %                 JA.noisyData_lie(i,i_targ) = 1; 
%         %             end
%         %         end
% 
%             end
%         end
%     end
    %% Calculate center of mass trajectory, store in JA
    if(saveJAStruct) % only worth execution time if saving data
        tic;
        dataLengthMax = max([size(JA.targAlign(1).jump(1).data,1), ...
                             size(JA.targAlign(2).jump(1).data,1), ...
                             size(JA.targAlign(3).jump(1).data,1)]);
        JA.CoMTraj = zeros(36,dataLengthMax,3);
        for i_targ = 1:numel(JA.targAlign)
            for i_jump = 1:numel(JA.targAlign(i_targ).jump)
                CoMTraj = getCoMTraj(JA,i_targ,i_jump,EKFCodePath);
                JA.CoMTraj(12*(i_targ-1)+i_jump,1:size(CoMTraj,1),:) = CoMTraj;
            end
        end
        disp('Time calculating CoM trajectories: ');
        toc
    end
    

    if(plotJointAngles)
        if(~use_shldrPrismModel)
            plotJoints = {[7,8],[10,15],[11,16],[12,17],[13,18],[20,27],[23,30],[25,32]};
%             plotJointsALL = {[7,8,9],[10,15],[11,16],[12,17],[13,18],[20,27],[21,28],[22,29],[23,30],[24,31],[25,32],[26,33]};
        else
            plotJoints = {[7,8],[11,17],[12,18],[13,19],[14,20],[22,29],[25,31],[27,34]};
        end
        plotJA_multi(JA,plotJoints,1);
    end


    if(saveJAStruct)
        mydir = pwd;
        idcs = strfind(mydir,'\');
        newdir = mydir(1:idcs(end)); 
        saveFilePath = [newdir 'results\JA_P' partNum '.mat'];

        save(saveFilePath,'JA');

        disp(['Participant ' partNum ' data struct saved'])
    end
    
    disp(['Participant ' partNum ' is complete.']);
    pause(1);
end


% flightMean(:,4) = mean(flightMean(:,1:3),2);
% flightVar(:,4) = mean(flightVar(:,1:3),2);
% flightVar_norm = flightVar./flightMean;
% [R,P] = corrcoef([partStability,flightMean(:,4),flightVar(:,4),flightVar_norm(:,4)]);


