function [param, traj_load, segmentInfo] = setup_main(filesToLoad, manSegLoadPath, currFilestack, mode, runSettings, cost_function_names, feature_load_names, extraSettings)
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

switch param.dataset
    case 'pr2'
%         filepathModel = 'project_pr2/robot.xml';
        filepathModel = 'robot.xml';
        mdl = rlCModel(filepathModel);
        mdl.forwardPosition();
        
%         bagFileLoad = convertCharsToStrings(filesToLoad.fullPathMat);
%         loadFeatureSet = process_data(bagFileLoad);
        loadMat = load(filesToLoad.fullPathMat);
        
        dt = 1/100;
        
%         loadFeatureSet.t = 0:dt:size(loadMat.POSITION, 2);
        loadFeatureSet.t = loadMat.TIME';
        loadFeatureSet.q = loadMat.POSITION';
        loadFeatureSet.dq = loadMat.VELOCITY';
        loadFeatureSet.ddq = loadMat.ACCELERATION';
        
      
       
        param.dofsFull = 1:7;
        param.dofsFromFull = 1:7;
        
        param.endEffectorName = '';        
        param.win_shift = external_win_shift;
        param.NbSample = size(loadFeatureSet.q, 1);          
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

param.model = mdl;

param.end_eff_frames = {'frame13'};
param.end_eff_framesInds = {1:3};

if 0
    vis = rlVisualizer('vis',640,480);
    vis.addModel(mdl);
    
    for i=1:numel(traj_load.t)
        mdl.position(param.dofsFromFull) = traj_load.q(:, i);
        mdl.velocity(param.dofsFromFull) = traj_load.dq(:, i);
        mdl.acceleration(param.dofsFromFull) = traj_load.ddq(:, i);
        
        mdl.forwardPosition();
        mdl.inverseDynamics();
        
        vis.update();
        pause(dt);
    end
end

param.win_shift = external_win_shift;
param = update_length(param, external_win_length, external_spline_length);

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


