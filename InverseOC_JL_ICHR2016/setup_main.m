function [param, traj_load, segmentInfo] = setup_main(filesToload, manSegLoadPath, currFilestack, mode, runSettings, cost_function_names, feature_load_names, extraSettings)
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
addOptional(p, 'doc_pivot_max_iter', 5000);

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

addOptional(p, 'optThreshold_doc', 1.745329251994e-9);
addOptional(p, 'optThreshold_ioc', 1.745329251994e-9);
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
external_optThreshold_doc = p.Results.optThreshold_doc;
external_optThreshold_ioc = p.Results.optThreshold_ioc;

param.corrThreshold = p.Results.correlation_threshold; % if the correlation is higher than this, then the higher entry is dropped

param.dataset = currFilestack.dataset;

param.jointChainStart = 'ankle';
param.datasetTag = runSettings.variableFactors.datasetTag;
param.outputString = runSettings.variableFactors.outputString;

param.doc_sim.knot_count = external_knots_sim;      % starting knot count (Adina: 11 knots to 100), previously using 4
param.doc_sim.maxiter = 5000;
param.doc_sim.TolFun = external_optThreshold_doc;
param.doc_sim.TolX = external_optThreshold_doc;
param.doc_sim.TolCon = external_optThreshold_doc;

param.ioc.knot_count = external_knots;    
param.ioc.maxiter = 5000;
param.ioc.TolFun = external_optThreshold_ioc;
param.ioc.TolX = external_optThreshold_ioc;
param.ioc.TolCon = external_optThreshold_ioc;

param.doc_pivot.knot_count = external_knots;     
param.doc_pivot.maxiter = external_doc_pivot_max_iter;
param.doc_pivot.TolFun = external_optThreshold_doc;
param.doc_pivot.TolX = external_optThreshold_doc;
param.doc_pivot.TolCon = external_optThreshold_doc; 


algorithmParam = structAlgorithmParam();
algorithmParam.ekfRun = [];
algorithmParam.linkDefinition = 'X00';       
algorithmParam.visualize = 0;
algorithmParam.endEffectorName = 'frame_rhand_end';

switch param.dataset
    case 'expressive_ioc'
        %filepathModel = '../../../PoseEstimation/kalmanfilter/ik_framework/instance_expressiveioc/model/ioc_v4_rightarm_fixedbase.xml';
        filepathModel = '../../../kalmanfilter/ik_framework/instance_expressiveioc/model/ioc_v4_rightarm_fixedbase.xml';
        %filepathModel = 'D:/aslab_gitlab/kalmanfilter/ik_framework/instance_expressiveioc/model/ioc_v4_rightarm_fixedbase.xml';

            % initialize the model instance and load marker data
            %     modelInstance = rlModelInstance_expressiveioc(currFileEntry.subjectNumber);
        modelInstance = rlModelInstance_expressiveioc_rightArm('');
        modelInstance.loadModelFromModelSpecsNoSensor(filepathModel, filesToload.fullPathMat);
        
%         dataInstance_trc = rlDataInstance_trc(modelInstance);
%         dataInstance_trc.loadData(filesToload.fullPathTrc, algorithmParam);
%         dataInstance_trc.dataProcessing(algorithmParam);

        % link the data instance into the model and create model to init XML model
        % linkage length and sensor attachments
%         [kinematicTransform, dynamicTransform, sensorTransform] = ...
%             modelInstance.makeModel(filepathModel, dataInstance_trc, algorithmParam);
        
        loadFeatureSet = rlFeatureSet_ioc();
        loadFeatureSet.loadDataFromFile(filesToload.fullPathMat);
        
        % resample data to 100 Hz
        loadFeatureSet.resample(1/100);
        loadFeatureSet.q = loadFeatureSet.q(:, 4:end); % remove the prismatic joint
        
        param.dofsFull = 1:7;
        param.dofsFromFull = 1:7;
        
        param.endEffectorName = algorithmParam.endEffectorName;

        dt = loadFeatureSet.dt;
        
        param.win_shift = external_win_shift;
        param.NbSample = size(loadFeatureSet.q, 1);
                
        if ~isempty(fieldnames(extraSettings)) && extraSettings.use
            pickFrame = round(extraSettings.pickTarget/2);
            placementFrame = round(extraSettings.placementTarget/2);
        end
end


% time params
param.ti_full = 0;
param.dt_full = dt;
        
% switch mode
%     case 'sim'        
%         param.tf_full = param.doc_sim_win_length - param.dt_full;
% 
%     case 'win'
%         loadFeatureSet.q = loadFeatureSet.q';
        param.tf_full = (size(loadFeatureSet.q, 1)-param.dt_full)*param.dt_full;
        
        if external_dataLoadLength == 0
            param.dataLoadLength = [];
        else
            param.dataLoadLength = external_dataLoadLength;
        end
% end

traj_load.t = 0:dt:param.tf_full;
fullDataAngles = loadFeatureSet.q;
traj_load.q =  filter_dualpassBW(fullDataAngles, 0.05, 0, 4)';
traj_load.dq = calcDeriv(traj_load.q, dt);
traj_load.ddq = calcDeriv(traj_load.dq, dt);
traj_load.dddq = calcDeriv(traj_load.ddq, dt);

segmentInfo.timeStart(1) = traj_load.t(1);
segmentInfo.timeEnd(1) = traj_load.t(end);

param.model = modelInstance.model;

param.end_eff_frames = {'frame_relbow_4', 'frame_rwrist_4', 'frame_rhand_end'};
param.end_eff_framesInds = {1:3, 4:6, 7:9};

massWeights = [param.model.bodies(2).m param.model.bodies(3).m param.model.bodies(4).m];
massWeights = massWeights / sum(massWeights);

param.weights.quantity_motion_dx = massWeights;
param.weights.weight_effort = massWeights;
param.weights.time_effort = massWeights;
param.weights.space_effort = massWeights;
param.weights.flow_effort = massWeights;
param.weights.shapeDir = massWeights;


if 0
    vis = rlVisualizer('vis',640,480);
    vis.addModel(modelInstance.model);
    
    for i=pickFrame:placementFrame%1:numel(traj_load.t)
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


if ~isempty(fieldnames(extraSettings)) && extraSettings.use
    % Get pick pose
    modelInstance.model.position(param.dofsFromFull) = traj_load.q(:, pickFrame);
    modelInstance.model.forwardPosition();
    pickPose = modelInstance.model.getFrameByName(algorithmParam.endEffectorName).t;

    % Get placement pose
    modelInstance.model.position(param.dofsFromFull) = traj_load.q(:, placementFrame);
    modelInstance.model.forwardPosition();
    placementPose = modelInstance.model.getFrameByName(algorithmParam.endEffectorName).t;

else
    
    strsplitStr = strsplit(filesToload.id, '_');
    subject = strsplitStr(1);
    
    switch subject{1}
        case 'Subject01'
            pickPose = [0.9405    0.3285   -0.0873    0.1286; 
                       -0.1124    0.0582   -0.9920    0.3439;
                       -0.3207    0.9427    0.0917   -0.2761; 
                             0         0         0    1.0000];
            placementPose = [0.8497    0.2090    0.4842    0.0965; 
                             0.4769    0.0872   -0.8746    0.3725;
                            -0.2250    0.9740   -0.0256   -0.1028; 
                                  0         0         0    1.0000];
        case 'Subject02'
            pickPose = [0.6467    0.2697   -0.7135    0.2186;
                       -0.7055   -0.1440   -0.6939    0.3974;
                       -0.2899    0.9521    0.0971   -0.3756;
                             0         0         0    1.0000];
           placementPose = [0.8758    0.4688    0.1150    0.0002;
                            0.1286    0.0031   -0.9917    0.4112;
                           -0.4652    0.8833   -0.0575   -0.3017;
                                 0         0         0    1.0000];
        case 'Subject03'
            pickPose = [0.7349    0.1300   -0.6656    0.3103;
                       -0.6484   -0.1529   -0.7458    0.3072;
                       -0.1988    0.9796   -0.0281   -0.4120;
                             0         0         0    1.0000];
            placementPose = [0.6273    0.7425    0.2352   -0.1200;
                             0.3798   -0.0280   -0.9246    0.3304;
                            -0.6799    0.6693   -0.2996   -0.3209;
                                  0         0         0    1.0000];
        case 'Subject04'
            pickPose = [0.9692    0.0888   -0.2298    0.1751;
                       -0.2463    0.3487   -0.9043    0.4004;
                       -0.0002    0.9330    0.3598   -0.3279;
                             0         0         0    1.0000];
            placementPose = [0.7780   -0.2995    0.5522   -0.0199;
                             0.6274    0.3247   -0.7078    0.3944;
                             0.0327    0.8971    0.4405   -0.2490;
                                  0         0         0    1.0000];
        case 'Subject05'
            pickPose = [0.7948    0.3007   -0.5272    0.0691;
                       -0.5238   -0.0990   -0.8461    0.3810;
                       -0.3066    0.9486    0.0788   -0.5646;
                             0         0         0    1.0000];
            placementPose = [0.6586    0.4822    0.5777   -0.1408;
                             0.4816    0.3198   -0.8160    0.3307;
                            -0.5782    0.8156   -0.0216   -0.4527;
                                  0         0         0   1.0000];
    end
%     % Default locations
%     param.x_anchor_1 = [-0.39075143  0.03671429 -0.46354857]'; 
%     param.x_anchor_2 = [-0.40905143 -0.04970571 -0.40449714]';
%     param.rot_anchor_1 = eye(3);
%     param.rot_anchor_2 = eye(3);
end


% Target locations
param.x_anchor_1 = pickPose(1:3, 4);
param.rot_anchor_1 = pickPose(1:3, 1:3);
param.x_anchor_2 = placementPose(1:3, 4);
param.rot_anchor_2 = placementPose(1:3, 1:3);

% Subwindow size for expressive parameters
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


