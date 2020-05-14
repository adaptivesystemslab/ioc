% Reads in previously-saved trc and ekf marker position data, replays
% visualization using ekf_state and lg_ekf_hatInvState (joint angles)

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
addpath(EKFCodePath);


dataFolder_JA = ['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data_2D\'];
dataFolder_DOC = ['C:\Users\kgwester\Documents\ResearchWork\ioc_project\jump_results\jump2d\test\'];

% if DOC ran cropped frame, list start and end frames
sf = 801;
ef = 1300;


% go two folder levels into saved DOC data
files = dir(dataFolder_DOC);
newFolderNum = find(contains({files.name},'2018'));
files = dir([dataFolder_DOC files(newFolderNum(1)).name '\']);
newFolderNum = find(contains({files.name},'P'));
jumpID = files(newFolderNum(1)).name
dataFolder_DOC = [files(newFolderNum(1)).folder '\' files(newFolderNum(1)).name '\'];
% find .mat file in folder
files = dir(dataFolder_DOC);
newFolderNum = find(contains({files.name},'.mat'));
dataFile_DOC = files(newFolderNum(1)).name

% load DOC data
load([dataFolder_DOC dataFile_DOC]);

% load original corresponding JA data struct
partNum = jumpID(2:3);
targNum = str2num(jumpID(6));
jumpNum = str2num(jumpID(9));

load([dataFolder_JA 'JA_P' partNum '_2D.mat']);




%Create original and DOC models
addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\scripts'));
[mdl_orig, ~] = createJumpModel_ioc_2D(JA, targNum, jumpNum, EKFCodePath, 0);
mdl_orig.forwardPosition;

[mdl_DOC, ~] = createJumpModel_ioc_2D(JA, targNum, jumpNum, EKFCodePath, 1);
% shift DOC model 1m to the left
mdl_DOC.transforms(1).t(1:3,4) = mdl_DOC.transforms(1).t(1:3,4) + [0,1,0]';
mdl_DOC.forwardPosition;



% Make trajectory files [DOFs x Frames]
q_orig = saveVar.feature_full.q;
q_DOC = [];
% for ind_windowCount = 1:numel(saveVar.q_recon_plot_array)
q_belowThreshold_combined = horzcat(saveVar.q_recon_plot_array{:});
t_belowThreshold_combined = horzcat(saveVar.t_recon_plot_array{:});

for ind_tCount = 1:length(saveVar.feature_full.t)
    currT = saveVar.feature_full.t(ind_tCount);     
    findTimeStep = find(t_belowThreshold_combined == currT);
    q_DOC(:, ind_tCount) = mean(q_belowThreshold_combined(:, findTimeStep), 2);
end


% Make Visualizer
vis = rlVisualizer('vis',640,960);
vis.addModel(mdl_orig);
vis.addModel(mdl_DOC);
vis.update();

for i= 1:size(q_orig,2)-1

    if(mod(i,1)==0)
        disp(['Frame: ' num2str(i)]);
            
        mdl_orig.position = q_orig(:,i);
        mdl_orig.forwardPosition;
        mdl_DOC.position = q_DOC(:,i);
        mdl_DOC.forwardPosition;
        
        vis.update();
    end

end

pause(1);

clear vis;