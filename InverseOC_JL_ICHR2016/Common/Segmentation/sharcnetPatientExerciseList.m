function [patientList, exerciseList] = sharcnetPatientExerciseList(dataSet, patientSet, movementSet)

switch dataSet
    case 'Healthy1'
        switch patientSet
            case 'One'
                for i = 1:20
                    patientList{i} = i;
                end
                
            case 'Group'
                patientList = {[2 3 4 6 8], [9 15 16 19 20], [5 7 10 11 12], [1 13 14 17 18]};
                
            case 'Group_healthy'
                patientList = {[2 3 4 6 8], [9 15 16 19 20], [5 7 10 11 12], [1 13 14 17 18]};
                
            case 'Group_test'
                patientList = {[1], [9 15 16 19 20], [5 7 10 11 12], [1 13 14 17 18]};
                
            case 'Group_training4'
                patientList = {[9 15 16 19 20]};
                
            case 'All'
                patientList = 1:20;
                
            otherwise
                patientList = {patientSet};
        end % trianingPatientProfile
        
        switch movementSet
            case 'All' % used for training
                exerciseStrings{1} = {'HFEO_SUP'};
                exerciseStrings{2} = {'KEFO_SIT'};
                exerciseStrings{3} = {'KHEF_SUP'};
                exerciseStrings{4} = {'SQUA_STD'};
                exerciseStrings{5} = {'STSO_SIT'};
                
                exerciseStrings = [exerciseStrings exerciseStrings];
                exerciseInstanceList = cell(1, 10);
                
                for i = 1:5
                    exerciseInstanceList{i} = 1;
                    exerciseInstanceList{i+5} = 2;
                end
                
            case 'All_OneCell' % used for testing
                exerciseStrings{1} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT', ...
                    'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT'};
                
                exerciseInstanceList{1} = [1 1 1 1 1 2 2 2 2 2];
                
            case 'Group' % spaced out to match with
                [exerciseStrings, exerciseInstanceList] = healthy1_group;
                
            case 'Group_healthy' % spaced out to match with healthy2 entries
                [exerciseStrings, exerciseInstanceList] = healthy1_group_healthy;
                
            otherwise
                exerciseStrings{1} = {movementSet};
                
                exerciseInstanceList{1} = zeros(size(exerciseStrings{1}));
        end % trainingmovementprofile
        
    case 'Healthy2'
        switch patientSet
            case 'One'
                for i = 1:10
                    patientList{i} = i;
                end
                
            case 'Group'
                patientList = {1:5, 6:10};
                
            case 'Group_healthy'
                patientList = {1:5, 1:5, 6:10, 6:10};
                
            case 'All'
                patientList = 1:10;
                
            case 'Group_training4'
                patientList = {[1:5]};
        end % patientset
        
        switch movementSet
            case 'All'
                exerciseStrings{1} = {'HAAO_STD'};
                exerciseStrings{2} = {'HAAO_SUP'};
                exerciseStrings{3} = {'HEFO_STD'};
                exerciseStrings{4} = {'HFEO_STD'};
                exerciseStrings{5} = {'KFEO_STD'};
                exerciseStrings{6} = {'KHEF_STD'};
                exerciseStrings{7} = {'LUNG_STD'};
                
                exerciseInstanceList = cell(1, 7);
                
                for i = 1:7
                    exerciseInstanceList{i} = 0;
                end
                
            case 'All_OneCell' % used for testing
                exerciseStrings{1} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};
                
                exerciseInstanceList{1} = [0 0 0 0 0 0 0];
                
            case 'Group'
                [exerciseStrings, exerciseInstanceList] = healthy2_group;
                
            case 'Group_healthy'
                [exerciseStrings, exerciseInstanceList] = healthy2_group_healthy;
        end % movementset
        
    case 'TRI1'
        switch patientSet
            case 'One'
                for i = 1:18
                    patientList{i} = i;
                end
                
            case 'Group_halfSeg'
                patientList = 6:10;
                
            case 'Group'
                patientList = {1:5, 6:10, 11:15, [1:3 16:18]};
                
            case 'Group_healthy'
                patientList = {1:5, 6:10, 11:15, [1:3 16:18]};
                
            case 'Group_training4'
                patientList = {[1:5]};
                
            case 'All'
                patientList = 1:18;
        end % patientset
        
        switch movementSet
            case 'All'
                exerciseStrings{1} = {'HFEO_SUP'};
                exerciseStrings{2} = {'KEFO_SIT'};
                exerciseStrings{3} = {'KHEF_SUP'};
                exerciseStrings{4} = {'SQUA_STD'};
                exerciseStrings{5} = {'STSO_SIT'};
                
                exerciseStrings{6} = {'KHEF_STD'};
                exerciseStrings{7} = {'HAAO_STD'};
                exerciseStrings{8} = {'HAAO_SUP'};
                exerciseStrings{9} = {'HEFO_STD'};
                exerciseStrings{10} = {'HFEO_STD'};
                exerciseStrings{11} = {'KFEO_STD'};
                exerciseStrings{12} = {'LUNG_STD'};
                
                for i = 1:13
                    exerciseInstanceList{i} = 0;
                end
                
            case 'All_OneCell' % used for testing
                exerciseStrings{1} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT', ...
                    'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};
                
                exerciseInstanceList{1} = zeros(1, 13);
                
            case 'KEFO_SIT'
                exerciseStrings{1} = {'KEFO_SIT'};
                exerciseInstanceList{1} = [0];
                
            case 'Group'
                [exerciseStrings, exerciseInstanceList] = tri_group;
                
            case 'Group_healthy'
                [exerciseStrings, exerciseInstanceList] = tri_group_healthy;
        end % movementset
        
    case 'StJoseph1'
        switch patientSet
            case 'One'
                for i = 1:26
                    patientList{i} = i;
                end
                
            case 'Group_training4'
                patientList = {[1:5]};
                
            case 'All'
                patientList = 1:26;
        end % patientset
        
    case 'TRO'
        switch patientSet
            case 'One'
                for i = 1:1
                    patientList{i} = i;
                end
        end % patientset
        
        switch movementSet
            case 'All'
                fullString = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
                    'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
                    'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
                    'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
                    'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'};
                
                for i = 1:length(fullString)
                    exerciseStrings{i} = fullString(i);
                end
                
                for i = 1:length(fullString)
                    exerciseInstanceList{i} = 0;
                end
                
            case 'All_OneCell' % used for testing
                exerciseStrings{1} = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
                    'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
                    'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
                    'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
                    'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'};
                
                exerciseInstanceList{1} = zeros(1, length(exerciseStrings{1}));
                
                
            case 'Group'
                [exerciseStrings, exerciseInstanceList] = tro_group;
        end
        
    case 'Myo'
        switch patientSet
            case 'One'
                for i = 1:10
                    patientList{i} = i;
                end
            case 'Group'
                patientList = {[1 3 5 7 9], [2 4 6 8 10]};
                
            case 'All'
                patientList = 1:10;
                
        end % patientset
        
        switch movementSet
            case 'All' % used for training
                exerciseStrings{1} = {'FIST_NON'};
                exerciseStrings{2} = {'FISP_NON'};
                exerciseStrings{3} = {'GUNM_NON'};
                exerciseStrings{4} = {'PONT_IND'};
                exerciseStrings{5} = {'PONT_INM'};
                exerciseStrings{6} = {'PDDM_INO'};
                exerciseStrings{7} = {'PDDM_OUT'};
                exerciseStrings{8} = {'SNAP_NON'};
                exerciseStrings{9} = {'THPK_NON'};
                
                %                 exerciseStrings = [exerciseStrings exerciseStrings];
                
                exerciseInstanceList = cell(1, 9);
                
                for i = 1:9
                    exerciseInstanceList{i} = 0;
                    %                     exerciseInstanceList{i+9} = [3 4 5];
                end
                
            case 'All_OneCell' % used for testing
                exerciseStrings{1} = {'FIST_NON','FISP_NON','GUNM_NON','PONT_IND',...
                    'PONT_INM','PDDM_INO','PDDM_OUT','SNAP_NON','THPK_NON','RAND_NON'};
                
                exerciseInstanceList{1} = [0 0 0 0 0 0 0 0 0 0];
                
            case 'All_OneCell_NoRand' % used for testing
                exerciseStrings{1} = {'FIST_NON','FISP_NON','GUNM_NON','PONT_IND',...
                    'PONT_INM','PDDM_INO','PDDM_OUT','SNAP_NON','THPK_NON'};
                
                exerciseInstanceList{1} = [0 0 0 0 0 0 0 0 0];
                
                
            case 'Rand' % used for testing
                exerciseStrings{1} = {'RAND_NON'};
                
                exerciseInstanceList{1} = [0];
                
            case 'Group' % spaced out to match with
                [exerciseStrings, exerciseInstanceList] = myo_group;
        end % trainingmovementprofile
        
    case 'Squats_TUAT_2011'
        switch patientSet
            case 'Training'
                patientList = {2:3};
                
            case 'Testing'
                patientList = 1:6;                
        end % patientset
        
        switch movementSet
            case {'All', 'SQUA_STD'} % used for training
                exerciseStrings{1} = {'SQUA_STD'};
                
                exerciseInstanceList{1} = 0;
        end % trainingmovementprofile
        
    case 'Squats_TUAT_2015'
        switch patientSet
            case 'Training'
                patientList = {2:3};
                
            case 'Testing'
                patientList = 1:6;        
            
            case 'Training_group'
                patientList{1} = [2 3 4 5 6 7 8];
                patientList{2} = [1 3 4 5 6 7 8];
                patientList{3} = [1 2 4 5 6 7 8];
                patientList{4} = [1 2 3 5 6 7 8];
                patientList{5} = [1 2 3 4 6 7 8];
                patientList{6} = [1 2 3 4 5 7 8];
                patientList{7} = [1 2 3 4 5 6 8];
                patientList{8} = [1 2 3 4 5 6 7];
                
            case 'Testing_group'
                patientList = 1:8;               
        end % patientset
        
        switch movementSet
            case {'All', 'SQUA_STD'} % used for training
                exerciseStrings{1} = {'SQUA_STD'};
                
                exerciseInstanceList{1} = 0;
        end % trainingmovementprofile
        
    case 'Taiso_UT_2009'
         switch patientSet
            case 'Training_Squat'
                patientList{1} = 1;
                sessionList{1} = {[20081209]};
                
                patientList{2} = 1;
                sessionList{2} = {[20090106]};
                
                patientList{3} = 1;
                sessionList{3} = {[20090203]};
                
                patientList{4} = 1;
                sessionList{4} = {[20090310]};
                
            case 'Testing_Squat'
                patientList = 1;  
                sessionList = {[20081209, 20090106, 20090203, 20090310]};
                
             case 'Training_Taiso'
                 patientList{1} = 1;
%                  sessionList{1} = {[20081209, 20090106, 20090203, 20090310]};
%                 sessio

             case 'Testing_Taiso'
                 patientList = 2:4;
         end
        
         switch movementSet
             case {'All', 'SQUA_STD'} % used for training
                 exerciseStrings{1} = {'SQUA_STD'};
                 
                 exerciseInstanceList{1} = 0;
                 
             case 'TAIS_STD_ONE'
                 exerciseStrings{1} = {'TAIS_STD_ONE'};
                 
                 exerciseInstanceList{1} = 0;
                 
             case 'TAIS_STD_TWO'
                 exerciseStrings{1} = {'TAIS_STD_TWO'};
                 
                 exerciseInstanceList{1} = 0;
         end % trainingmovementprofile
        
    case 'Doppel'
        switch patientSet
            case 'Training'
                patientList = {1};
                
            case 'Testing'
                patientList = 1;
        end % patientset
        
        switch movementSet
            case 'Training'
                exerciseStrings{1} = {'SQUA_STD_NON'};
                
                exerciseInstanceList{1} = 1;
            case 'Testing'
                exerciseStrings{1} = {'SQUA_STD_NON'};
                
                exerciseInstanceList{1} = 2;
                
        end % trainingmovementprofile
end % trainingDataset

% master list of group activities (if have not been covered already)
switch movementSet
    case 'Group_top5'
        [exerciseStrings, exerciseInstanceList] = group_top5;
        
    case 'Group_top5_group1'
        [exerciseStrings, exerciseInstanceList] = group_top5_group1;
        
    case 'Group_top5_group2'
        [exerciseStrings, exerciseInstanceList] = group_top5_group2;
        
    case 'Group_top5_group12'
        [exerciseStrings, exerciseInstanceList] = group_top5_group12;
        
    case 'Group_top5_group3'
        [exerciseStrings, exerciseInstanceList] = group_top5_group3;
        
    case 'Group_second5'
        [exerciseStrings, exerciseInstanceList] = group_second5;
end

if ~exist('sessionList', 'var')
    sessionList = cell(size(exerciseInstanceList));
end

exerciseList.sessionList = sessionList;
exerciseList.exerciseString = exerciseStrings;
exerciseList.exerciseInstance = exerciseInstanceList;
exerciseList.exerciseCount = length(exerciseStrings);

end

function [exerciseList, exerciseInstanceList] = healthy1_group
exerciseList{1} = {'KEFO_SIT', 'STSO_SIT'};
exerciseList{2} = {'SQUA_STD'};
exerciseList{3} = {'HFEO_SUP', 'KHEF_SUP'};
exerciseList{4} = {'KEFO_SIT', 'HFEO_SUP', 'KHEF_SUP'};
exerciseList{5} = {'STSO_SIT', 'SQUA_STD'};
exerciseList{6} = {'KEFO_SIT'};
exerciseList{7} = {'HFEO_SUP'};
exerciseList{8} = {'STSO_SIT', 'KHEF_SUP', 'SQUA_STD'};
exerciseList{9} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT'};

exerciseList{10} = {'KEFO_SIT', 'STSO_SIT'};
exerciseList{11} = {'SQUA_STD'};
exerciseList{12} = {'HFEO_SUP', 'KHEF_SUP'};
exerciseList{13} = {'KEFO_SIT', 'HFEO_SUP', 'KHEF_SUP'};
exerciseList{14} = {'STSO_SIT', 'SQUA_STD'};
exerciseList{15} = {'KEFO_SIT'};
exerciseList{16} = {'HFEO_SUP'};
exerciseList{17} = {'STSO_SIT', 'KHEF_SUP', 'SQUA_STD'};
exerciseList{18} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT'};

for i = 1:9
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*1;
    exerciseInstanceList{i+9} = ones(size(exerciseList{i}))*2;
end

end

function [exerciseList, exerciseInstanceList] = healthy1_group_healthy
exerciseList{1} = {'KEFO_SIT', 'STSO_SIT'}; % ip - sit
exerciseList{2} = {'SQUA_STD'};             % ip - std
exerciseList{3} = {'HFEO_SUP', 'KHEF_SUP'}; % ip - sup
exerciseList{4} = {'KEFO_SIT', 'HFEO_SUP', 'KHEF_SUP'}; % kc - ha
exerciseList{5} = {'STSO_SIT', 'SQUA_STD'};             % kc - ah
exerciseList{6} = {'KEFO_SIT'};                         % aj - knee
exerciseList{7} = {'HFEO_SUP'};                         % aj - hip
exerciseList{8} = {'STSO_SIT', 'KHEF_SUP', 'SQUA_STD'}; % aj - both
exerciseList{9} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT'};  % all - healthy1
exerciseList{10} = {'None'};                                                     % all - healthy2
exerciseList{11} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT'}; % all - all

exerciseList = [exerciseList exerciseList];

for i = 1:11
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*1;
    exerciseInstanceList{i+11} = ones(size(exerciseList{i}))*2;
end

end

function [exerciseList, exerciseInstanceList] = healthy2_group
exerciseList{1} = {'HAAO_STD', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};
exerciseList{2} = {'HAAO_SUP'};
exerciseList{3} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};
exerciseList{4} = {'KFEO_STD'};
exerciseList{5} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD'};
exerciseList{6} = {'KHEF_STD', 'LUNG_STD'};
exerciseList{7} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};

for i = 1:7
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = healthy2_group_healthy
exerciseList{1} = {'None'};                                                                 % ip - sit
exerciseList{2} = {'HAAO_STD', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; % ip - std
exerciseList{3} = {'HAAO_SUP'};                                                             % ip - sup
exerciseList{4} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; % kc - ha
exerciseList{5} = {'None'};                                                                             % kc - ah
exerciseList{6} = {'KFEO_STD'};                                        % aj - knee
exerciseList{7} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD'};    % aj - hip
exerciseList{8} = {'KHEF_STD', 'LUNG_STD'};                            % aj - both
exerciseList{9} = {'None'};                                                                               % all - healthy1
exerciseList{10} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};  % all - healthy2
exerciseList{11} = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'};  % all - all

exerciseList = [exerciseList exerciseList];

for i = 1:11
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
    exerciseInstanceList{i+11} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = tri_group
% basically, a combined version of healthy1_group and healthy2_group
exerciseList{1} = {'KEFO_SIT', 'STSO_SIT'}; % sit
exerciseList{2} = {'SQUA_STD', 'HAAO_STD', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; % stand
exerciseList{3} = {'HFEO_SUP', 'KHEF_SUP', 'HAAO_SUP'}; % supine
exerciseList{4} = {'KEFO_SIT', 'HFEO_SUP', 'KHEF_SUP', 'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; %h2a
exerciseList{5} = {'STSO_SIT', 'SQUA_STD'}; % a2h
exerciseList{6} = {'KEFO_SIT', 'KFEO_STD'}; % knee
exerciseList{7} = {'STSO_SIT', 'KHEF_SUP', 'SQUA_STD', 'KHEF_STD', 'LUNG_STD'}; % hip and knee
exerciseList{8} = {'HFEO_SUP', 'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD'}; % hip
exerciseList{9} = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT', 'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; % all

for i = 1:9
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList]  = tri_group_healthy
[exercise_healthy1, exerciseInstanceList_healthy1] = healthy1_group_healthy;
[exercise_healthy2, exerciseInstanceList_healthy2] = healthy2_group_healthy;

for i = 1:length(exercise_healthy1)
    exerciseList{i} = [exercise_healthy1(i) exercise_healthy2(i)];
    exerciseInstanceList{i} = [exerciseInstanceList_healthy1(i) exerciseInstanceList_healthy2(i)];
end
end

function [exerciseList, exerciseInstanceList] = group_top5
%     specStruct.exercise = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'}; % type in the movement names here. please be consistent with string length, ie {'KEFO_SUP', 'KHEF_SIT} or {'KEFO_SUP_NON1'}

exerciseList{1} = {'KHEF_SUP'};
exerciseList{2} = {'KEFO_SIT'};
exerciseList{3} = {'HEFO_STD'};
exerciseList{4} = {'KHEF_STD'};
exerciseList{5} = {'HAAO_STD'};

for i = 1:5
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = group_top5_group1
%     specStruct.exercise = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'}; % type in the movement names here. please be consistent with string length, ie {'KEFO_SUP', 'KHEF_SIT} or {'KEFO_SUP_NON1'}

exerciseList{1} = {'HEFO_STD', 'KHEF_STD', 'HAAO_STD'};

for i = 1
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = group_top5_group2
%     specStruct.exercise = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'}; % type in the movement names here. please be consistent with string length, ie {'KEFO_SUP', 'KHEF_SIT} or {'KEFO_SUP_NON1'}

exerciseList{1} = {'KHEF_SUP', 'KHEF_STD'};

for i = 1
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = group_top5_group12
%     specStruct.exercise = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'}; % type in the movement names here. please be consistent with string length, ie {'KEFO_SUP', 'KHEF_SIT} or {'KEFO_SUP_NON1'}

exerciseList{1} = {'HEFO_STD', 'KHEF_STD', 'HAAO_STD'};
exerciseList{2} = {'KHEF_SUP', 'KHEF_STD'};

for i = 1:2
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = group_top5_group3
%     specStruct.exercise = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'}; % type in the movement names here. please be consistent with string length, ie {'KEFO_SUP', 'KHEF_SIT} or {'KEFO_SUP_NON1'}

% % % exerciseList{1} = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'};

exerciseList{1} = {'HAAO_STD', 'HAAL_STD', 'HEFO_STD', 'HFEO_STD', 'KEFO_SIT', 'KFEO_SIT', 'KHEF_STD', 'KHEF_SUP'};

for i = 1
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = group_second5
%     specStruct.exercise = {'KHEF_SUP', 'KEFO_SIT', 'HEFO_STD', 'KHEF_STD', 'HAAO_STD'}; % type in the movement names here. please be consistent with string length, ie {'KEFO_SUP', 'KHEF_SIT} or {'KEFO_SUP_NON1'}
%  KFEO_SUP 'SQUA_STD', 'HFEO_STD', 'KEFO_SUP', 'KFEO_STD'

% % % exerciseList{1} = {'KFEO_SUP', 'SQUA_STD', 'HFEO_STD', 'KEFO_SUP', 'KFEO_STD'};

exerciseList{1} = {'HAAO_SUP', 'HFEO_SUP', 'KFEO_STD', 'KEFO_SUP', 'KFEO_SUP', 'LUNG_STD', 'SQUA_STD'};

for i = 1
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = tro_group

% group by...
% 1. sets of motions  % not included: LER RRL
% 2. active joint angles
% 3. upper/lower body

exerciseList{1} = {'BAD', 'BAU'};
exerciseList{2} = {'BAL180', 'BAR180'};
exerciseList{3} = {'BAL90', 'BAR90'};
exerciseList{4} = {'LAL180', 'LAR180'};
exerciseList{5} = {'LAL90', 'LAR90'};
exerciseList{6} = {'LKAL', 'LKAR', 'LKE', 'LKR'};
exerciseList{7} = {'LPAL', 'LPAR', 'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR'};
exerciseList{8} = {'MLL', 'MLR', 'MRL', 'MRR'};
exerciseList{9} = {'RAL180', 'RAR180'} ;
exerciseList{10} = {'RAL90', 'RAR90'};
exerciseList{11} = {'RKAL', 'RKAR', 'RKE', 'RKR'};
exerciseList{12} = {'RPAL', 'RPAR', 'RPE', 'RPR'};
exerciseList{13} = {'SQD', 'SQU'};
exerciseList{14} = {'WLL', 'WLR', 'WRL', 'WRR'};

exerciseList{15} = {'LAL180', 'LAL90', 'LAR180', 'LAR90', 'LKAL', 'LKAR', 'LPAL', 'LPAR', 'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR'}; % left upper
exerciseList{16} = {'RAL180', 'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RPAL', 'RPAR', 'RPE', 'RPR'}; % right upper
exerciseList{17} = {'BAL180', 'BAL90', 'BAR180', 'BAR90'}; % all upper
exerciseList{18} = {'LKE', 'LKR'}; % left lower
exerciseList{19} = {'RKE', 'RKR'}; % right lower
exerciseList{20} = {'BAD', 'BAU', 'MLL', 'MLR', 'MRL', 'MRR', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'}; % all body

exerciseList{21} = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
    'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
    'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
    'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
    'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'};


for i = 1:21
    exerciseInstanceList{i} = ones(size(exerciseList{i}))*0;
end
end

function [exerciseList, exerciseInstanceList] = myo_group
exerciseList{1} = {'FIST_NON', 'GUNM_NON', 'PONT_IND', 'PONT_INM', 'SNAP_NON'}; % motion - full hand
exerciseList{2} = {'FISP_NON', 'THPK_NON'}; % motion - subtle
exerciseList{3} = {'PDDM_INO', 'PDDM_OUT'}; % motion - wrist
exerciseList{4} = {'FIST_NON', 'GUNM_NON', 'PONT_IND', 'PONT_INM', 'SNAP_NON', 'THPK_NON'}; % fingers - thumb
exerciseList{5} = {'FIST_NON', 'GUNM_NON', 'PONT_IND', 'PONT_INM', 'SNAP_NON'}; % fingers - index
exerciseList{6} = {'FIST_NON', 'GUNM_NON', 'PONT_INM', 'SNAP_NON'}; % fingers - middle
exerciseList{7} = {'FISP_NON', 'PDDM_INO', 'PDDM_OUT'}; % fingers - none
exerciseList{8} = {'FIST_NON', 'FISP_NON', 'GUNM_NON', 'PONT_IND', 'PONT_INM', 'PDDM_INO', 'PDDM_OUT', 'SNAP_NON', 'THPK_NON'}; % all

% exerciseList = [exerciseList exerciseList exerciseList exerciseList exerciseList];

upperLim = 8;
for i = 1:upperLim
    exerciseInstanceList{i+0*upperLim} = ones(size(exerciseList{i}))*0;
    %     exerciseInstanceList{i+1*upperLim} = ones(size(exerciseList{i}))*2;
    %     exerciseInstanceList{i+2*upperLim} = ones(size(exerciseList{i}))*3;
    %     exerciseInstanceList{i+3*upperLim} = ones(size(exerciseList{i}))*4;
    %     exerciseInstanceList{i+4*upperLim} = ones(size(exerciseList{i}))*5;
end

end