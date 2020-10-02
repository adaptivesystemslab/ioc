function [param, traj_load, segmentInfo] = setup_main(filesToload, manSegLoadPath, currFilestack, mode, runSettings, cost_function_names, feature_load_names)
% a common setup file for main_direct and main_inverse

% parse the input runSettings.varaibleFactors
p = inputParser;
p.KeepUnmatched = true;

addOptional(p, 'doc_sim_win_length', 2); % length of data generated from DOC [s]
addOptional(p, 'knots_sim', 4); % length of data generated from DOC [s]

addOptional(p, 'win_length', 61); % the window over the observed data [sample]
addOptional(p, 'win_shift', 20); % the slide of the above window [sample]
addOptional(p, 'spline_length', 61); % the window length of the spline itself. if mismatching from win_length, the data is interpolated [sample]
addOptional(p, 'spline_order', 5);
addOptional(p, 'knots', 10); % the knots used in IOC/DOC. regardless the knot count, all constraints are added as knots
addOptional(p, 'correlation_threshold', 1);

addOptional(p, 'dofsToUse', 7:45);
addOptional(p, 'pivotMethod', 'resnorm');
addOptional(p, 'pivotSelection', 'all');
addOptional(p, 'dataLoadLength', 0);

addOptional(p, 'half_phase_construction', 'no');
addOptional(p, 'doc_sim_sp_dq_endcond', 'none');
addOptional(p, 'doc_sim_cf_dq_constraint', 'constraint_dxddx');
addOptional(p, 'doc_pivot_max_iter', 500);

addOptional(p, 'ioc_gradient_method', 'numerical');
addOptional(p, 'ioc_cf_q_constraint', 'startmid1endonknot'); % should the constraints location be variable in the reconstruction? 3 points (fixed), matching constraint (match). 
addOptional(p, 'ioc_cf_q_constraint_y', 'previous'); 
addOptional(p, 'ioc_cf_dq_constraint', 'constraint_dxddx'); 
addOptional(p, 'ioc_cf_dq_constraint_y', 'previous_dxddx'); 
addOptional(p, 'ioc_knot_locations', 'variable');
addOptional(p, 'ioc_add_cf_as_knots', 'no'); 
addOptional(p, 'ioc_add_knot_as_cf_constraints', 'none'); % if 1, then all knots are also added as constraints
addOptional(p, 'ioc_sp_dq_endcond', 'none'); 

addOptional(p, 'doc_rmse_cf_q_constraint', 'startmid1endonknot'); % should the constraints location be variable in the reconstruction? 3 points (fixed), matching constraint (match). 
addOptional(p, 'doc_rmse_cf_q_constraint_y', 'previous'); 
addOptional(p, 'doc_rmse_cf_dq_constraint', 'constraint_dxddx'); 
addOptional(p, 'doc_rmse_cf_dq_constraint_y', 'previous_dxddx'); 

addOptional(p, 'optThreshold', 1.745329251994e-9);
addOptional(p, 'normalizationMethod_feature', 'none');
addOptional(p, 'normalizationMethod_cf', 'rang');

addOptional(p, 'splineType', 'piecewise_cds');
addOptional(p, 'splineSlave', '5th_order_poly');
addOptional(p, 'windowMode', 'wide');
addOptional(p, 'segment_only_windows', 'none'); % if yes, then win_length/spline_length/win_shift is ignored and adjusted according to some manual segment

addOptional(p, 'simDistance', 100);
addOptional(p, 'dynamicSampleRate', 0);

parse(p, runSettings.variableFactors);
        
param = p.Results;

external_win_length = p.Results.win_length; 
external_win_shift = p.Results.win_shift;
external_spline_length = p.Results.spline_length;
external_knots = p.Results.knots;
external_knots_sim = p.Results.knots_sim;
external_doc_pivot_max_iter = p.Results.doc_pivot_max_iter;
external_dataLoadLength = p.Results.dataLoadLength;
external_optThreshold = p.Results.optThreshold;

param.corrThreshold = p.Results.correlation_threshold; % if the correlation is higher than this, then the higher entry is dropped

param.dataset = currFilestack.dataset;

param.jointChainStart = 'ankle';
param.datasetTag = runSettings.variableFactors.datasetTag;
param.outputString = runSettings.variableFactors.outputString;

tolFun = 1.745329251994e-9; % corresponds to deg2rad(1e-6)
tolX =   1.745329251994e-9;
tolCon = 1.745329251994e-9;

param.doc_sim.knot_count = external_knots_sim;      % starting knot count (Adina: 11 knots to 100), previously using 4
param.doc_sim.maxiter = 1000;
param.doc_sim.TolFun = tolFun;
param.doc_sim.TolX = tolX;
param.doc_sim.TolCon = tolCon;

param.ioc.knot_count = external_knots;    
param.ioc.maxiter = 500;
param.ioc.TolFun = tolFun;
param.ioc.TolX = tolX;
param.ioc.TolCon = tolCon;

tolFun = external_optThreshold; 
tolX = external_optThreshold;
tolCon = external_optThreshold;

param.doc_pivot.knot_count = external_knots;     
param.doc_pivot.maxiter = external_doc_pivot_max_iter;
param.doc_pivot.TolFun = tolFun;
param.doc_pivot.TolX = tolX;
param.doc_pivot.TolCon = tolCon; 


algorithmParam = structAlgorithmParam();
% algorithmParam = structAlgorithmParam;
algorithmParam.visualize = 0;
algorithmParam.endEffectorName = 'frame_t1c7_0';

switch param.dataset
    case 'iit_2017'
        filepathModel = fullfile('.', 'Models', 'iit_v8.xml');
        
        loadFeatureSet = rlFeatureSet_ioc();
        loadFeatureSet.loadDataFromFile(filesToload{1});
        
        modelInstance = rlModelInstance_iit(currFilestack.subjectNumber);
        modelInstance.loadMarkerData(filesToload{2});
        modelInstance.makeModel(filepathModel);
        
        param.dofsFull = loadFeatureSet.inds_pos;
        param.dofsFromFull = loadFeatureSet.inds_model_pos;
        
        param.endEffectorName = 'frame_t1c7_0';
        
    case 'iit_2017_2d'
        filepathModel = fullfile('.', 'Models', 'iit_v8_2d.xml');
        
        loadFeatureSet = rlFeatureSet_ioc();
        loadFeatureSet.loadDataFromFile(filesToload{1});
        
        modelInstance = rlModelInstance_iit_2d(currFilestack.subjectNumber);
        modelInstance.loadMarkerData(filesToload{2});
        modelInstance.makeModel(filepathModel);
        
        param.dofsFull = loadFeatureSet.inds_pos;
        param.dofsFromFull = loadFeatureSet.inds_model_pos;
        
        %         loadFeatureSet.q(:, [27:29 35:37 43:45] - 6) = 0;
        %         indsVias = [918 1048 1189];
        
        param.endEffectorName = 'frame_t1c7_0';
        
    case 'squats_urfi_2011'
        loadFeatureSet = load(filesToload{1}); % q, dq, ddq, Markers
        load3 = load(filesToload{2}); % body_height, body_weight, FP, COP

        dt = 0.01;
        fullDataAngles = loadFeatureSet.q;
        loadFeatureSet.q = dualpassButterworth(fullDataAngles', 0.05, 0, 4);
        
        dofsToUse = 2:4;
        
        q = loadFeatureSet.q(:, dofsToUse);
        
        Markers = loadFeatureSet.Markers;
        param.body_weight = load3.body_weight;
        param.body_height = load3.body_height;
        
        segmentInfoLoad = arsLoader(manSegLoadPath);
        segmentInfo.use = segmentInfoLoad.use;
        segmentInfo.segmentCount = 1:length(segmentInfoLoad.use);
        segmentInfo.timeStart = segmentInfoLoad.timeStart / 1000;  % convert frames to seconds
        segmentInfo.timeEnd = segmentInfoLoad.timeEnd / 1000;
        
        param.L1=abs(Markers.RANK(1,1)); % RANK, frame 1, X
        param.L2=abs(Markers.RANK(1,2)); % RANK, frame 1, Y
        param.L3=norm(Markers.RANK(1,1:2)-Markers.RKNE(1,1:2)); % RANK-RKNE, frame 1, X&Y
        param.L4=norm(Markers.RKNE(1,1:2)-Markers.RHJC(1,1:2));
        param.L5=norm(Markers.RHJC(1,1:2)-Markers.Mid_Pelvis(1,1:2)); % hip to pelvis (from ASIS/PSIS)
        param.L6=0; 
        param.L7=0;
        param.L8=0;    
        param.L9=0;
        param.L10=0;
        param.L11=0;
        
        fatigue_level = [];
        
        param.dofsFull = 1:7;
        param.dofsFromFull = dofsToUse;
        
        param.win_shift = external_win_shift;
        param.NbSample = size(q, 2);
        
%         param.ti_full = 0;
%         param.dt_full = 1/100;
%         param.tf_full = (size(q, 2)-param.dt_full)*param.dt_full;
        loadFeatureSet.time = (0:0.01:size(q, 2))*dt;
        
        param.endEffectorName = 'frame5';
        
        filepathModel = fullfile('.', 'Models', '7dof_squat.xml');
        modelInstance.model = rlCModel(filepathModel);
        modelInstance.model.forwardPosition();
        param = setup_dyn(param);
        modelInstance.model = updateModelInfo(modelInstance.model, param);
        
    case 'healthy1ankle'
        filepathModel = fullfile('.', 'Models', 'healthy1ankle_3joints_111.xml');
        
        subjectParam.exerciseName = 'SQUA_STD_SLO1';
        subjectParam.fileId = '';
        subjectParam.subjectNumber = 1;
        subjectParam.sessionNumber = 1;
        
        subjectParam.filePath = fullfile('D:', 'aslab', 'data', 'Lowerbody_healthy1_2011-11', ...
            ['Subject' num2str(subjectParam.subjectNumber)], ['Session' num2str(subjectParam.sessionNumber)], subjectParam.exerciseName);
        
        pathToKinematicData = fullfile('D:', 'aslab', 'data', 'IOCDatasets', 'healthy1_ankle_2d', 'data',  ...
            [subjectParam.exerciseName '_' 'Subj' num2str(subjectParam.subjectNumber) '_' 'Sess' num2str(subjectParam.sessionNumber) '.mat']);
    
        loadFeatureSet = rlFeatureSet_ioc();
        loadFeatureSet.loadDataFromFile(pathToKinematicData);
%         loadFeatureSet.cropLoadedData(1:3300);
%         loadFeatureSet.q(:, [27:29 35:37 43:45] - 6) = 0;
%         indsVias = [1062 1276 1663];
%         
% segmentInfo.timeStart(1) = traj_load.t(indsVias(1));
% segmentInfo.timeEnd(1) = traj_load.t(indsVias(3));

        param.endEffectorName = 'frame_shoulder0';
        
        modelInstance = rlModelInstance_healthy1ankle(subjectParam.subjectNumber);
        
        param.dofsFull = 1:3;
        param.dofsFromFull = 1:3;
        param.body_weight=60;
        param.body_height=1.7;
        
    case {'jump2D','jump2d'}
        %% Load planar jumping data
        cropFrames = 0; % crop calibration motion of each jump recording
        
        load(filesToload{1}); % Loads JA Struct: q, dq, ddq, part demographic info
%         load3 = load(filesToload{2}); % body_height, body_weight, FP, COP
        targNum = str2num(filesToload{2}(1));
        jumpNum = str2num(filesToload{2}(3:end));
        
        dt = 0.005;
        numJoints = numel(JA.jointNames);
        fullDataAngles = JA.targ(targNum).jump(jumpNum).data(1+cropFrames:end,:);
%         loadFeatureSet.q = dualpassButterworth(fullDataAngles', 0.05, 0, 4); % NEED TRANSPOSE??? 
        loadFeatureSet.q = fullDataAngles';
        
        dofsToUse = 1:numJoints;
        
        q = fullDataAngles(:, dofsToUse)';
        
%         Markers = loadFeatureSet.Markers;
        param.body_weight = (JA.partData.weight)/2.2; % kg conversion
        param.body_height = JA.partData.height;
        
        param.dofsFull = 1:numJoints;
        param.dofsFromFull = dofsToUse;
        
        param.win_shift = external_win_shift;
        param.NbSample = size(q, 2);
        
%         param.ti_full = 0;
%         param.dt_full = 1/100;
%         param.tf_full = (size(q, 2)-param.dt_full)*param.dt_full;
        loadFeatureSet.time = (0:size(q, 2))*dt;
        
        param.endEffectorName = 'rtoe0';
        
        % Create 2D jumping model
        addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\scripts'));
%         filepathModel = 'C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\scripts\JumpModel_Eul_inertia_2D.xml';
        EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
        [jump_2D_model, initPos] = createJumpModel_ioc_2D(JA, targNum, jumpNum, EKFCodePath);
        modelInstance.model = jump_2D_model;
        modelInstance.model.position = initPos;
        modelInstance.model.velocity(:) = 0;
        modelInstance.model.acceleration(:) = 0;
        modelInstance.model.forwardPosition;
        modelInstance.model.forwardVelocity;
%         param = setup_dyn(param);
%         modelInstance.model = updateModelInfo(modelInstance.model, param);
        
        
        segmentInfo.use = 1;
        
        
        % Set additional jump-specific parameters
        param.jump.takeoffFrame = JA.TOFrame(jumpNum,targNum);
        param.jump.landFrame = JA.LandFrame(jumpNum,targNum);
        param.jump.locationLand = JA.locationLand(12*(targNum-1) + jumpNum);
        param.jump.grade = JA.jumpGrades(12*(targNum-1) + jumpNum);
        param.jump.modelLinks = JA.modelLinks;
        param.jump.world2base = squeeze(JA.world2base((12*(targNum-1) + jumpNum),:,:));
        param.jump.bad2d = JA.bad2D(jumpNum,targNum);
        
        
        
end


% time params
param.ti_full = 0;
param.dt_full = dt;
        
% switch mode
%     case 'sim'        
%         param.tf_full = param.doc_sim_win_length - param.dt_full;
% 
%     case 'win'
        loadFeatureSet.q = loadFeatureSet.q';
        param.tf_full = (size(loadFeatureSet.q, 1)-param.dt_full)*param.dt_full;
        
        if external_dataLoadLength == 0
            param.dataLoadLength = [];
        else
            param.dataLoadLength = external_dataLoadLength;
        end
% end

traj_load.t = 0:dt:param.tf_full;
fullDataAngles = loadFeatureSet.q(:, param.dofsFromFull);
traj_load.q =  dualpassButterworth(fullDataAngles, 0.05, 0, 4)';
traj_load.dq = calcDeriv(traj_load.q, dt);
traj_load.ddq = calcDeriv(traj_load.dq, dt);
traj_load.dddq = calcDeriv(traj_load.ddq, dt);

segmentInfo.timeStart(1) = traj_load.t(1);
segmentInfo.timeEnd(1) = traj_load.t(end);

% Set segmentInfo.time parameters so jump flight section is seen as
% segment in plot_sharcnet
segmentInfo.timeStart = (param.jump.takeoffFrame - cropFrames)*dt;  % convert frames to seconds
segmentInfo.timeEnd = (param.jump.landFrame - cropFrames)*dt;

param.model = modelInstance.model;

if 0
    vis = rlVisualizer('vis',640,480);
    vis.addModel(modelInstance.model);
    
    for i=1:numel(traj_load.t)
        modelInstance.model.position(param.dofsFromFull) = traj_load.q(:, i);
        modelInstance.model.velocity(param.dofsFromFull) = traj_load.dq(:, i);
        modelInstance.model.acceleration(param.dofsFromFull) = traj_load.ddq(:, i);
        
        modelInstance.model.forwardPosition();
        modelInstance.model.inverseDynamics();
        
        vis.update();
        pause(dt);
    end
end

param.win_shift = external_win_shift;
param = update_length(param, external_win_length, external_spline_length);

param = setupNormalizationValues(param, traj_load, segmentInfo, param.normalizationMethod_feature, param.normalizationMethod_cf);

% % constraint points. the x (time) and y (angle) positions for these locks
const_x(1) = 1; % THIS IS INDEX OF t, NOT TIME
const_x(2) = ceil(length(param.t_spline)/2); 
const_x(3) = length(param.t_spline);

const_y(:, 1) = traj_load.q(:, 1) ;
const_y(:, 2) = traj_load.q(:, param.simDistance);
const_y(:, 3) = traj_load.q(:, 1) ;

% const_y(:, 1) = traj_load.q(:, indsVias(1));
% const_y(:, 2) = traj_load.q(:, indsVias(2));
% const_y(:, 3) = traj_load.q(:, indsVias(3));

param.const_x = const_x;
param.const_y = const_y;
param.const_dx = [];
param.const_dy = [];
param.const_ddx = [];
param.const_ddy = [];
param.dof_count = size(traj_load.q, 1);

% param = update_intermed_ind(param, param.doc_pivot.knot_count);

param.removeHighlyCorrelatedJ = 1;

% gradient params
param.h = 1e-6;
param.diffType = 'central';
param.recoverThreshold_c = 1e-12; % if the recovered weights is less than than this, it is zeroed
param.recoverThreshold_j = 1e-12; % if the recovered contribution ratio is less than this, it is zeroed

param.cost_function_names = cost_function_names;
param.feature_load_names = feature_load_names;


