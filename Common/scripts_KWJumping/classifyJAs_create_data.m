% load JA data from multiple participants, train LSTM on joint angle trajectories 
clear;


partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};


X_all = zeros(35,500,21*36); % [DoF x maxNumberOfFrames x allJumps]
Info_all = zeros(21*36,5); % [allJumps x (partNum+targNum+jumpNum+jumpGradeNumber+toeDistToTarg]
dataLength = zeros(21,3); % cropped, aligned length of sets of 12 jump trajectories

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
                
                X_all(:,1:dataLength(i_part, i_targ),36*(i_part-1)+12*(i_targ-1)+i_jump) = ...
                    JA.targAlign(i_targ).jump(i_jump).data(startFrame:endFrame,:)';
                
                jumpGradeLetter = JA.jumpGrades{12*(i_targ-1)+i_jump};
%                 jumpGradeNumber = find(contains(jumpGradeTypes,jumpGradeLetter)==1);
                switch jumpGradeLetter{1}
                    case 'B'
                        jumpGradeNumber = 1;
                    case 'SB'
                        jumpGradeNumber = 2;
                    case 'P'
                        jumpGradeNumber = 3;
                    case 'P*'
                        jumpGradeNumber = 4;
                    case 'SF'
                        jumpGradeNumber = 5;
                    case 'F'
                        jumpGradeNumber = 6;
                end
                
                toeDistToTarg = 0;
                Info_all(36*(i_part-1)+12*(i_targ-1)+i_jump,:) = ...
                    [i_part+1, i_targ, i_jump, jumpGradeNumber, toeDistToTarg];
                
            end
        end
    end
    
end


% remove extra zeros in X_all
maxDataLength = max(max(dataLength));
X_all = X_all(:,1:maxDataLength,:);


% convert X_all to cell array of observations (required for 2018a LSTM
% network)
X_all_tmp = [];
for i = 1:size(X_all,3)
    X_all_tmp{i} = X_all(:,:,i);
end
X_all = X_all_tmp;

% save data
saveFileName = 'JA_NN_data.mat';
save(saveFileName,'X_all','Info_all');
disp('JA_NN_data.mat saved');

