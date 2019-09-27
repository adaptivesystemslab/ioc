% Reads in previously-saved joint angle data and makes plots, compares
% joint angles of one participant over consecutive jumps (see "learning")

% Choose dataset and parameters

plotJointAngles = 0;
fixOutliers = 1;
saveJAStruct = 1;


partNum = '02';
targetDist = {'55','70','85'};

JA = []; %"joint angles all"
JA.partNum = partNum;
JA.missingData = zeros(12,3);
JA.noisyData_eul = zeros(12,3);
JA.noisyData_lie = zeros(12,3);
JA.shiftFlipData_eul = zeros(36,33);
JA.shiftFlipData_lie = zeros(36,33);
numJoints = 0;


%% Read in and set participant data and jump grading
[age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
JA.partData.age = age;
JA.partData.gender = gender;
JA.partData.height = height/100; %in [m]
JA.partData.weight = weight;
JA.partData.S2sc = S2sc/100; %in [m]
[targetLengths,jumpGrades] = getJumpGrading(partNum);
JA.targetLengths = targetLengths;
JA.jumpGrades = jumpGrades;


%% Read in joint angle data for all jumps for this participant
for i_targ = 1:3
    for i_set = 1:2
        for i_jump = 1:6
            currJump = 6*(i_set-1) + i_jump;
            
            mydir = pwd;
            idcs = strfind(mydir,'\');
            newdir = mydir(1:idcs(end)-1); 
            saveFilePath = [newdir '\results\RESULTS_JA_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
            if(exist(saveFilePath,'file')==2)
                load(saveFilePath);
                
                JA.targ_eul(i_targ).jump(currJump).data = jointAngles_eul;
                JA.targ_lie(i_targ).jump(currJump).data = jointAngles_lie;
                JA.dataLengths(currJump,i_targ) = size(jointAngles_eul,1);
                JA.startFrame(currJump,i_targ) = startFrame;
                JA.TOFrame(currJump,i_targ) = footTOFrame;
                JA.LandFrame(currJump,i_targ) = footLandFrame;
                JA.endFrame(currJump,i_targ) = endFrame;
                
                if(numJoints==0)
                    JA.jointNames = jointAngles_names;
                    numJoints = size(jointAngles_eul,2);
                end
                    
            else %file not available, populate JAall with zeros
                JA.missingData(6*(i_set-1) + i_jump,i_targ) = 1;
                
                if(numJoints==0)    numJoints = 33;     end %first file missing, assume 33 joints
                JA.targ_eul(i_targ).jump(currJump).data = zeros(2500,numJoints);
                JA.targ_lie(i_targ).jump(currJump).data = zeros(2500,numJoints);
                JA.dataLengths(currJump,i_targ) = 2500;
                JA.startFrame(currJump,i_targ) = 1300;
                JA.TOFrame(currJump,i_targ) = 1500;
                JA.LandFrame(currJump,i_targ) = 1600;
                JA.endFrame(currJump,i_targ) = 1800;
            end
        end
    end
end



%% Align all joint angle data to take-off frame, crop frames near start and end of each recording
JA.targAlign_eul = JA.targ_eul;
JA.targAlign_lie = JA.targ_lie;
dataLengthsAlign = JA.dataLengths;

for i_targ = 1:numel(JA.targAlign_eul)
    % crop frames off front to align "footTOFrame" for all 12 jumps
    [minTOFrame,idx_TOFrame] = min(JA.TOFrame(:,i_targ));
    JA.TOAlignJump(i_targ) = idx_TOFrame;
    
    for currJump = 1:numel(JA.targAlign_eul(i_targ).jump)
        if(currJump~=idx_TOFrame)
            cropFrames = JA.TOFrame(currJump,i_targ) - minTOFrame;
            dataLengthsAlign(currJump,i_targ) = dataLengthsAlign(currJump,i_targ) - cropFrames;
            
            JA.targAlign_eul(i_targ).jump(currJump).data = ...
                JA.targAlign_eul(i_targ).jump(currJump).data(cropFrames+1:end,:);
            JA.targAlign_lie(i_targ).jump(currJump).data = ...
                JA.targAlign_lie(i_targ).jump(currJump).data(cropFrames+1:end,:);
        end
    end
    
    % crop frames off back to make all 12 jumps the same length
    [minDataLength,idx_dataLength] = min(dataLengthsAlign(:,i_targ));
    for currJump = 1:12
        if(currJump~=idx_dataLength)
            JA.targAlign_eul(i_targ).jump(currJump).data = ...
                JA.targAlign_eul(i_targ).jump(currJump).data(1:minDataLength,:);
            JA.targAlign_lie(i_targ).jump(currJump).data = ...
                JA.targAlign_lie(i_targ).jump(currJump).data(1:minDataLength,:);
        end
    end
end


% vertically shift/flip data, find and remove outliers remove noisy data
for i_targ = 1:numel(JA.targAlign_eul)
    idx_TOFrame = JA.TOAlignJump(i_targ);
    % get all jump data for specific joint
    for j = 1:numJoints
        sig_eul = zeros(size(JA.targAlign_eul(i_targ).jump(1).data,1),numel(JA.targAlign_eul(i_targ).jump));
        sig_lie = zeros(size(JA.targAlign_lie(i_targ).jump(1).data,1),numel(JA.targAlign_lie(i_targ).jump));
        for i = 1:size(sig_eul,2)
            sig_eul(:,i) = JA.targAlign_eul(i_targ).jump(i).data(:,j);
            sig_lie(:,i) = JA.targAlign_lie(i_targ).jump(i).data(:,j);
        end
        
        % vertically shift/flip JA signals, remove outlier signals, and remove noisy signals
        if(fixOutliers) 
            alignWindow = [JA.startFrame(idx_TOFrame,i_targ),JA.endFrame(idx_TOFrame,i_targ)];
            [~,shiftSig_eul, shift2Sig_eul, flipSig_eul, shiftFlipSig_eul, shiftFlip2Sig_eul, OL_smooth_eul] ...
                = cleanAndAlignData(sig_eul,alignWindow,1);
            [~,shiftSig_lie, shift2Sig_lie, flipSig_lie, shiftFlipSig_lie, shiftFlip2Sig_lie, OL_smooth_lie] ...
                = cleanAndAlignData(sig_lie,alignWindow,1);
        end
        
        % fix signal shifting/outlier/noise removal for EUL data
        for i = 1:numel(JA.targAlign_eul(i_targ).jump)
            tmpData = JA.targAlign_eul(i_targ).jump(i).data(:,j);
            if(shift2Sig_eul(i)~=0)
                JA.targAlign_eul(i_targ).jump(i).data(:,j) = tmpData + 2*pi*shift2Sig_eul(i);
                JA.shiftFlipData_eul(12*(i_targ-1) + i,j) = 1;
            elseif(shiftSig_eul(i)~=0)
                JA.targAlign_eul(i_targ).jump(i).data(:,j) = tmpData + pi*shiftSig_eul(i);
                JA.shiftFlipData_eul(12*(i_targ-1) + i,j) = 1;
            elseif(shiftFlip2Sig_eul(i)~=0)
                JA.targAlign_eul(i_targ).jump(i).data(:,j) = -(tmpData + 2*pi*shiftFlip2Sig_eul(i));
                JA.shiftFlipData_eul(12*(i_targ-1) + i,j) = 1;
            elseif(shiftFlipSig_eul(i)~=0)
                JA.targAlign_eul(i_targ).jump(i).data(:,j) = -(tmpData + pi*shiftFlipSig_eul(i));
                JA.shiftFlipData_eul(12*(i_targ-1) + i,j) = 1;
            elseif(flipSig_eul(i)~=0)
                JA.targAlign_eul(i_targ).jump(i).data(:,j) = -tmpData;
                JA.shiftFlipData_eul(12*(i_targ-1) + i,j) = 1;
            end
            % NOTE: ordered signal mods above from least to most likely to have false positives
            
            if(OL_smooth_eul(i)~=0)
%                 JA.targAlign_eul(i_targ).jump(i).data(:,j) = zeros(size(tmpData));
                JA.noisyData_eul(i,i_targ) = 1; 
            end
        end
        
        % fix signal shifting/outlier/noise removal for LIE data
        for i = 1:numel(JA.targAlign_lie(i_targ).jump)
            tmpData = JA.targAlign_lie(i_targ).jump(i).data(:,j);
            if(shift2Sig_lie(i)~=0)
                JA.targAlign_lie(i_targ).jump(i).data(:,j) = tmpData + 2*pi*shift2Sig_lie(i);
                JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
            elseif(shiftSig_lie(i)~=0)
                JA.targAlign_lie(i_targ).jump(i).data(:,j) = tmpData + pi*shiftSig_lie(i);
                JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
            elseif(shiftFlip2Sig_lie(i)~=0)
                JA.targAlign_lie(i_targ).jump(i).data(:,j) = -(tmpData + 2*pi*shiftFlip2Sig_lie(i));
                JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
            elseif(shiftFlipSig_lie(i)~=0)
                JA.targAlign_lie(i_targ).jump(i).data(:,j) = -(tmpData + pi*shiftFlipSig_lie(i));
                JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
            elseif(flipSig_lie(i)~=0)
                JA.targAlign_lie(i_targ).jump(i).data(:,j) = -tmpData;
                JA.shiftFlipData_lie(12*(i_targ-1) + i,j) = 1;
            end
            % NOTE: ordered signal mods above from least to most likely to have false positives
            
            if(OL_smooth_lie(i)~=0)
%                 JA.targAlign_lie(i_targ).jump(i).data(:,j) = zeros(size(tmpData));
                JA.noisyData_lie(i,i_targ) = 1; 
            end
        end
        
    end
end


if(plotJointAngles)
    plotJoints = {[7,8,9],[10,15],[11,16],[12,17],[13,18],[20,27],[23,30],[25,32],[26,33]};
    plotJointsALL = {[7,8,9],[10,15],[11,16],[12,17],[13,18],[20,27],[21,28],[22,29],[23,30],[24,31],[25,32],[26,33]};
    plotJA_multi(JA,plotJoints,1,1);
end


if(saveJAStruct)
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)); 
    saveFilePath = [newdir 'results\JA_P' partNum '.mat'];
    
    save(saveFilePath,'JA');

    disp(['Participant ' partNum ' data struct saved'])
end

