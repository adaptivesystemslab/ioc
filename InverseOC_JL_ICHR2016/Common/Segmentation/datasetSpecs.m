function specs = datasetSpecs(dataset)
    % this function holds the dataset name, the filepath root, and the full
    % exercise set usable for each dataset. in some of these cases, the
    % folder path and the dataset name is not set to be the same. in those
    % cases, both names are listed 
    
    % listing both the dataset name and its filepath name, in case they're
    % different
    switch lower(dataset)
        case {'healthy1', 'lowerbody_healthy1_2011-11'} % lowerbodyrehab_UW_2011
            exerciseListFull = {'KEFO_SIT_SLO', 'KHEF_SUP_SLO', 'SQUA_STD_SLO', 'STSO_SIT_SLO', 'HFEO_SUP_SLO'}; % full list
            datasetName = 'healthy1';
            dataPathSuffix = 'Lowerbody_healthy1_2011';
            
        case {'healthy1ankle'} % lowerbodyrehab_UW_2011
            exerciseListFull = {'STSO_SIT_SLO', 'HFEO_SUP_SLO'}; % full list
            datasetName = 'healthy1ankle';
            dataPathSuffix = 'Lowerbody_healthy1_2011';
            
        case {'healthy1hip'} % lowerbodyrehab_UW_2011
            exerciseListFull = {'KEFO_SIT_SLO', 'KHEF_SUP_SLO', 'SQUA_STD_SLO'}; % full list
            datasetName = 'healthy1hip';
            dataPathSuffix = 'Lowerbody_healthy1_2011';            
            
        case {'healthy2', 'lowerbody_healthy2_2013-07'} % lowerbodyrehab_UW_2013
            exerciseListFull = {'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; % full list
            datasetName = 'healthy2';
            dataPathSuffix = 'Lowerbody_healthy2_2013';
            
        case {'stjoseph1', 'lowerbody_stjoseph1_2013-02'} % lowerbodyrehab_sjhc_2013
            exerciseListFull = {'KEFO_SIT'};
            datasetName = 'stjoseph1';
            dataPathSuffix = 'Lowerbody_StJoseph1_2013';
            
        case {'tri1', 'lowerbody_tri1_2012-10'} % lowerbodyrehab_tri_2012
            exerciseListFull = {'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT', 'HFEO_SUP', ...
                'HAAO_STD', 'HAAO_SUP', 'HEFO_STD', 'HFEO_STD', 'KFEO_STD', 'KHEF_STD', 'LUNG_STD'}; % full list
            datasetName = 'tri1';
            dataPathSuffix = 'Lowerbody_TRI1_2012';
            
        case {'kulic2009_tro'} % fullbody_kulic_2009_tro
            exerciseListFull = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
                'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
                'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
                'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
                'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'}; % full list
            datasetName = 'kulic2009_tro';
            dataPathSuffix = 'kulic2009_tro';
            
        case {'myo', 'thalmic_2014-12'}
            exerciseListFull = {'FIST_NON', 'FISP_NON', 'GUNM_NON', 'PONT_IND', ...
                'PONT_INM', 'PDDM_INO', 'PDDM_OUT', 'SNAP_NON', 'THPK_NON', 'RAND_NON'}; % full list
            datasetName = 'myo';
            dataPathSuffix = 'Thalmic_2014-12';  
            
        case {'squats_tuat_2011', 'squats_tuat_2011-06'} % squats_tuat_2011
            exerciseListFull = {'SQUA_STD'};
            datasetName = 'squats_tuat_2011';
            dataPathSuffix = 'Squats_TUAT_2011';
            
        case {'squats_urfi_2011', 'squats_tuat_2015', 'squats_tuat_2015-12'} % squats_tuat_2011
            exerciseListFull = {'SQUA_STD'};
            datasetName = 'squats_urfi_2011';
            dataPathSuffix = 'Squats_URFI_2011';
            
        case {'doppel'} % squats_tuat_2011
            exerciseListFull = {'SQUA_STD'};
            datasetName = 'doppel';
            dataPathSuffix = 'Doppel';
            
        case {'taiso_ut_2009'} % taiso_ut_2009
            exerciseListFull = {'TAIS_STD_ONE', 'TAIS_STD_TWO'};
%             exerciseListFull = {...
%                 'MoveBody3Times', 'ArmCircleInnerOuter', 'SideBendingRightLeft', 'FrontBending', ...
%                 'WaistRotationRightLeft', 'LegsArms', 'TouchFootLeftRight', 'BigCirclesRightLeft', 'Jump'}; % TAISOU_1
            datasetName = 'taiso_ut_2009';
            dataPathSuffix = 'Taiso_UT_2009';
            
        case 'iit_2017'
            exerciseListFull = {'SQUA_STD'};
            datasetName = 'iit_2017';
            dataPathSuffix = fullfile('FullBody_IIT_2017', '2018_09_25_5', 'matEkfIk');
            
        case 'iit_2017_2d'
            exerciseListFull = {'SQUA_STD'};
            datasetName = 'iit_2017_2d';
            dataPathSuffix = 'Fullbody_IIT_2017_IK';
            
        case 'healthy1ik'
            exerciseListFull = {'SQUA_STD_SLO'};
            datasetName = 'healthy1ik';
            dataPathSuffix = 'Lowerbody_healthy1_2011_IK';
    end

    specs.exerciseListFull = exerciseListFull;
    specs.datasetName = datasetName;
    specs.dataPathSuffix = dataPathSuffix;
end