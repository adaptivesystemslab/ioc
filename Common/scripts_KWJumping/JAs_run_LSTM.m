% load JA data from multiple participants, train LSTM on joint angle trajectories 
clear;

% NN LSTM parameters
inputSize = 35; % number of joints
classType = 'jumpGrade'; % jumpGrade, toeDistToTarg, custom
numClasses = 3; % too short, on target, too long

numHiddenUnits = 100; % hiddeon units of LSTM, EXPERIMENT WITH THIS

partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'04'};


% target_percents = zeros(numel(partsToRun),2);


X_all = zeros(35,500,21*36); % [DoF x maxNumberOfFrames x allJumps]
Y_all = zeros(21*36,5); % [allJumps x (partNum+targNum+jumpNum+jumpGradeNumber+toeDistToTarg]
dataLength = zeros(21*3); % cropped, aligned length of sets of 12 jump trajectories

jumpGradeTypes = {'B','SB','P','P*','SF','F'}; %this is numbering order of "jumpGradeNumber" in Y_all

for i_part = 1:numel(partsToRun)
    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partsToRun{i_part} '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' i_part '.mat NOT FOUND']);
    else
        load(loadFilePath); % JA data struct
        
        
        
        for i_targ = 1:3
            endFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ) + 10; % extra 0.05 seconds to make sure feet have landed
            startFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ) - 300;
            dataLength(i_part, i_targ) = endFrame - startFrame + 1;
            
            for i_jump = 1:12
                
                X_all(:,1:dataLength,36*(i_part-1)+12*(i_targ-1)+i_jump) = ...
                    JA.targAlign(i_targ).jump(i_jump).data(startFrame:endFrame,:)';
                
                jumpGradeLetter = JA.jumpGrades{12*(i_targ-1)+i_jump};
                jumpGradeNumber = find(contains(jumpGradeTypes,jumpGradeLetter)==1);
                toeDistToTarg = 0;
                Y_all(36*(i_part-1)+12*(i_targ-1)+i_jump,:) = ...
                    [i_part+1, i_targ, i_jump, jumpGradeNumber, toeDistToTarg];
                
            end
        end
    end
    
end







