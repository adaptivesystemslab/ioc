%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
clear;

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);

saveModelLinkLengths = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'08'};


for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};


    data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));
    
    
    modelLinks_rec = [];
    modelLinks_ave = [];
    world2base_rec = zeros(36,4,4);
    
    
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
                    [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_v2(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);

                    mdl_eul.forwardPosition;
                    %Shift the Lie group model by 1 meter in x and y to visualize better
    %                 mdl_lie.transforms(1).t(1:2,4) = mdl_eul.transforms(1).t(1:2,4) + [1 1]';
    %                 mdl_lie.forwardPosition();

    
                    
                    currJump = 12*(i_targ-1) + 6*(i_set-1) + i_jump;
                    modelLinks_rec.jump(currJump) = modelLinks;
                    world2base_rec(currJump,:,:) = mdl_eul.transforms(1).t;
                    
                    modelLinks_ave.spine(currJump) = modelLinks.spine;
                    modelLinks_ave.spine2shldr_L(currJump,1:3) = modelLinks.spine2shldr_L;
                    modelLinks_ave.spine2shldr_R(currJump,1:3) = modelLinks.spine2shldr_R;
                    modelLinks_ave.base2hip_L(currJump,1:3) = modelLinks.base2hip_L;
                    modelLinks_ave.base2hip_R(currJump,1:3) = modelLinks.base2hip_R;
                    modelLinks_ave.uparm_L(currJump) = modelLinks.uparm_L;
                    modelLinks_ave.forearm_L(currJump) = modelLinks.forearm_L;
                    modelLinks_ave.uparm_R(currJump) = modelLinks.uparm_R;
                    modelLinks_ave.forearm_R(currJump) = modelLinks.forearm_R;
                    modelLinks_ave.thigh_L(currJump) = modelLinks.thigh_L;
                    modelLinks_ave.shin_L(currJump) = modelLinks.shin_L;
                    modelLinks_ave.foot_L(currJump,1:3) = modelLinks.foot_L;
                    modelLinks_ave.thigh_R(currJump) = modelLinks.thigh_R;
                    modelLinks_ave.shin_R(currJump) = modelLinks.shin_R;
                    modelLinks_ave.foot_R(currJump,1:3) = modelLinks.foot_R;
                    
                    
                    
%                     if(saveModelLinkLengths)
%                         mydir = pwd;
%                         idcs = strfind(mydir,'\');
%                         newdir = mydir(1:idcs(end)); 
%                         saveFilePath = [newdir 'results\ModelLinkLengths_JA_P' partNum '.mat'];
%                         
%                         world2base = mdl_eul.transforms(1).t;
%                         
%                         save(saveFilePath,'jointAngles_names','jointAngles_eul',...
%                                 'footTOFrame','footLandFrame','startFrame','endFrame',...
%                             'modelLinks','world2base','locationStart','locationLand');
%     %                     save(saveFilePath,'jointAngles_names','jointAngles_eul','jointAngles_lie',...
%     %                             'footTOFrame','footLandFrame','startFrame','endFrame');
% 
%                         disp([dataName ' JA data saved'])
%                     end


                    disp([dataName ' model complete'])
                end
    % END MAIN NESTED FOR LOOPS
            end
        end
    end

    % finished this participant
    
    
    modelLinkFields = fieldnames(modelLinks_ave);
    for i = 1:numel(modelLinkFields)
        modelLinks_ave.(modelLinkFields{i}) = mean(modelLinks_ave.(modelLinkFields{i}));
    end
    
    
    
    if(saveModelLinkLengths)
        mydir = pwd;
        idcs = strfind(mydir,'\');
        newdir = mydir(1:idcs(end)); 
        saveFilePath = [newdir 'results\Model_Info_JA_P' partNum '.mat'];

        save(saveFilePath,'modelLinks_rec','modelLinks_ave','world2base_rec');

        disp(['P' partNum ' JA model info saved'])
    end
end

