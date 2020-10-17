function [specStruct, variableFactors] = script3_testingSpec
% set up base settings
setting_costfunction = '8';
setting_dataset = 'pr2';
setting_iocvalidation = 'resnorm';
setting_outputstring = 'test';
setting_normalization = 'none'; % rang rset none

sourceDataMatFolder = 'D:/data/data_pr2';
resultsDataFolder = 'D:/data/results_pr2';



variableFactors.outputPathSuffix = [setting_normalization '_' setting_costfunction '_' 'cfSumSq'];

docsim_gen_settings_20knots; % 10knots 20knots
variableFactors.knots = 20; % 4 7 10 12 // 10 20
cconst_array = [1 1 1 1 1 1];

switch setting_costfunction
    case '7'
        cost_function_names{1} = 'half_joint_task';
        cost_function_names{2} = 'joint_length';
        cost_function_names{3} = 'manip_rot';
        cost_function_names{4} = 'manip_trans';
        cost_function_names{5} = 'manipulability'; %'';
        cost_function_names{6} = 'orientation_length'; %
        cost_function_names{7} = 'task_length'; %
        
        ccost_array{1} = [1 0 0 0 0 0 0 0];
        
    case '8'
        cost_function_names{1} = 'half_joint_task';
        cost_function_names{2} = 'joint_length';
        cost_function_names{3} = 'joint_limit';
        cost_function_names{4} = 'manip_rot';
        cost_function_names{5} = 'manip_trans';
        cost_function_names{6} = 'manipulability'; %'';
        cost_function_names{7} = 'orientation_length'; %
        cost_function_names{8} = 'task_length'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0];    
end

variableFactors.outputString = setting_outputstring;
variableFactors.win_length = 41; % 61 121 181 -- 101 201
variableFactors.win_shift = 10;
variableFactors.ioc_cf_q_constraint = 'startmid1endonknotbelow'; %  (fixed, variable, edge, startmidend, startmidmidend) where to set the consts? 3 points (fixed), global matching constraint (variable).

variableFactors.generateSimOnly = 0;
variableFactors.spline_length = variableFactors.win_length;

variableFactors.correlation_threshold = 1; % [0:1]
variableFactors.splineType = 'piecewise_cds'; % bformspline cubicspline polyfit ppformspline slm splinefit_lundgren piecewise_cds
variableFactors.splineSlave = '5th_order_poly'; % if the master spline is piecewise, what should the secondary one be?

variableFactors.half_phase_construction = 'no'; % yes no double
variableFactors.doc_sim_sp_dq_endcond = 'none'; % (constraint_dx, two_dx, two_dxddx, none) for traj gen, where should the ECC spline constraints be?
variableFactors.doc_sim_cf_dq_constraint = 'constraint_dxddx'; % (two_dx, constraint_dx, none) where should sp constraints be?

variableFactors.ioc_gradient_method = 'numerical'; % numerical, symbolic
variableFactors.ioc_cf_q_constraint_y = 'previous'; % (fixed, previous) where to set the consts? 3 points (fixed), global matching constraint (variable).
variableFactors.ioc_cf_dq_constraint = 'constraint_dxddx'; % (edge_dx variable edge_dxddx) where to set the consts? 3 points (fixed), global matching constraint (variable).
variableFactors.ioc_cf_dq_constraint_y = 'previous_dxddx'; % (fixed_dx, previous_dx, zero_dx) where to set the consts? 3 points (fixed), global matching constraint (variable).
variableFactors.ioc_knot_locations = 'variable'; % where to set the knots in the IOC? (balanced, variable)
variableFactors.ioc_add_knot_as_cf_constraints = 'none'; % (no, yes) if yes, then all normal knots are also added as constraints
variableFactors.ioc_sp_dq_endcond = 'none'; % () which points are the sp constrains
variableFactors.ioc_add_cf_as_knots = 'no';

variableFactors.doc_rmse_cf_q_constraint = variableFactors.ioc_cf_q_constraint; % (variableEdge fixed, variable, edge) where to set the consts? 3 points (fixed), global matching constraint (variable).
variableFactors.doc_rmse_cf_q_constraint_y = variableFactors.ioc_cf_q_constraint_y; % (previous fixed) where to set the consts? 3 points (fixed), global matching constraint (variable).
variableFactors.doc_rmse_cf_dq_constraint = variableFactors.ioc_cf_dq_constraint; % (variableEdge_dxddx edge_dx variable edge_dxddx) where to set the consts? 3 points (fixed), global matching constraint (variable).
variableFactors.doc_rmse_cf_dq_constraint_y = variableFactors.ioc_cf_dq_constraint_y; % (previous_dxddx fixed_dx, previous_dx, zero_dx) where to set the consts? 3 points (fixed), global matching constraint (variable).

variableFactors.segment_only_windows = 'none'; % (none, constraint, seg, mergeseg) if yes, then win_length/spline_length/win_shift is ignored and adjusted according to some manual segment
variableFactors.windowMode = 'wide'; % wide narrow

specStruct.ccost_array = ccost_array;
specStruct.cconst_array = cconst_array;

specStruct.dataset = setting_dataset;
specStruct.patient = [1];
specStruct.exerciseAcceptPrefixString = {''};
specStruct.exerciseAcceptPrefixInstance = 0;
runMode = 'win';

variableFactors.datasetTag = specStruct.dataset;

variableFactors.normalizationMethod_feature = 'none';
variableFactors.normalizationMethod_cf = setting_normalization;
   
variableFactors.pivotMethod = setting_iocvalidation;
variableFactors.pivotSelection = 'all'; 
variableFactors.optThreshold = deg2rad(1e-6); % corresponds to deg2rad(1e-6)

variableFactors.simDistance = 100;
variableFactors.dynamicSampleRate = 0;

variableFactors.batchWindowBreakdownFactor = 3000; % (3000) break down the batch runs to increments of this factor, \pm 2 seconds on both sides
variableFactors.batchWindowBreakdownFactorOverlap = 600;

main_batch(specStruct, runMode, cost_function_names, setting_outputstring, variableFactors, sourceDataMatFolder, resultsDataFolder);
