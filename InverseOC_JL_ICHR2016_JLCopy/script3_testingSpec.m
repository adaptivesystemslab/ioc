function [specStruct, variableFactors] = script3_testingSpec
% set up base settings
setting_costfunction = '8';
setting_dataset = 'expressive_ioc';
setting_iocvalidation = 'resnorm';
setting_outputstring = 'test';

% setting_pathToRawData = fullfile('D:', 'aslab', 'data');
% setting_pathToOutput = fullfile('D:', 'results');

setting_pathToRawData = '../../expressiveiocData/';
setting_pathToOutput = '../../expressiveiocData/results';

docsim_gen_settings_20knots; % 10knots 20knots
variableFactors.knots = 20; % 4 7 10 12 // 10 20
cconst_array = [1 1 1 1 1 1];

switch setting_costfunction
    case '3'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'ddx';
        cost_function_names{3} = 'tau';
        
        ccost_array{1} = [1 0 0];
        ccost_array{2} = [0 1 0];
        ccost_array{3} = [0 0 1];
        ccost_array{4} = [1 1 1];
       
    case '8'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'dddx';
        cost_function_names{4} = 'tau';
        cost_function_names{5} = 'dtau';
        cost_function_names{6} = 'ddtau'; %'';
        cost_function_names{7} = 'dq_tau'; %
        cost_function_names{8} = 'kinetic_en_sqrt'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0];    
     
    case '12_base'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'ddx';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'x_cartCurv';
        cost_function_names{6} = 'tau';
        cost_function_names{7} = 'dtau';
        cost_function_names{8} = 'ddtau'; %'';
        cost_function_names{9} = 'dq_tau'; %
        %cost_function_names{10} = 'x_anchor_1'; % Distance to pick location
        %cost_function_names{11} = 'rot_anchor_1'; % Distance to pick hand orientation
        cost_function_names{10} = 'x_anchor_2'; % Distance to placement location 
        cost_function_names{11} = 'rot_anchor_2'; % Distance to placement hand orientation
        cost_function_names{12} = 'kinetic_en_sqrt'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0 0 0 0 0];               
        
    
    case '14_base'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'ddx';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'x_cartCurv';
        cost_function_names{6} = 'tau';
        cost_function_names{7} = 'dtau';
        cost_function_names{8} = 'ddtau'; %'';
        cost_function_names{9} = 'dq_tau'; %
        cost_function_names{10} = 'x_anchor_1'; % Distance to pick location
        cost_function_names{11} = 'rot_anchor_1'; % Distance to pick hand orientation
        cost_function_names{12} = 'x_anchor_2'; % Distance to placement location 
        cost_function_names{13} = 'rot_anchor_2'; % Distance to placement hand orientation
        cost_function_names{14} = 'kinetic_en_sqrt'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0 0 0 0 0 0 0];               
        
    case '14'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'dddx';
        cost_function_names{4} = 'tau';
        cost_function_names{5} = 'dtau';
        cost_function_names{6} = 'ddtau'; %'';
        cost_function_names{7} = 'dq_tau'; %
        cost_function_names{8} = 'kinetic_en_sqrt'; % 
        cost_function_names{9} = 'quantity_motion_dx'; % 
        cost_function_names{10} = 'volume_bounding_box'; % 
        cost_function_names{11} = 'weight_effort'; % 
        cost_function_names{12} = 'time_effort'; % 
        cost_function_names{13} = 'space_effort'; % 
        cost_function_names{14} = 'flow_effort'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
        
     case '15'
        cost_function_names{1} = 'tau';
        cost_function_names{2} = 'dq_tau'; %
        cost_function_names{3} = 'kinetic_en_sqrt'; % 
        cost_function_names{4} = 'quantity_motion_dx'; % 
        cost_function_names{5} = 'volume_bounding_box'; % 
        cost_function_names{6} = 'weight_effort'; % 
        cost_function_names{7} = 'time_effort'; % 
        cost_function_names{8} = 'space_effort'; %     
        cost_function_names{9} = 'x_anchor_1'; % 
        cost_function_names{10} = 'x_anchor_2'; % 
        cost_function_names{11} = 'com'; % 
        cost_function_names{12} = 'x_displace'; %
        cost_function_names{13} = 'dx_displace'; %
        cost_function_names{14} = 'cartCurv'; % 
        cost_function_names{15} = 'dcom'; % 
        
        ccost_array{1} = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
        
     case '17'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'dddx';
        cost_function_names{4} = 'tau';
        cost_function_names{5} = 'dtau';
        cost_function_names{6} = 'ddtau'; %'';
        cost_function_names{7} = 'dq_tau'; %
        cost_function_names{8} = 'kinetic_en_sqrt'; % 
        cost_function_names{9} = 'quantity_motion_dx'; % 
        cost_function_names{10} = 'volume_bounding_box'; % 
        cost_function_names{11} = 'weight_effort'; % 
        cost_function_names{12} = 'time_effort'; % 
        cost_function_names{13} = 'space_effort'; % 
        cost_function_names{14} = 'flow_effort'; % 
     
        cost_function_names{15} = 'x_anchor_1'; % 
        cost_function_names{16} = 'x_anchor_2'; % 
        cost_function_names{17} = 'com'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];      
     
     case '17_final'
         
        cost_function_names{1} = 'ddx';
        cost_function_names{2} = 'tau';
        cost_function_names{3} = 'dq_tau'; %
        cost_function_names{4} = 'kinetic_en_sqrt'; % 
        cost_function_names{5} = 'quantity_motion_dx'; % 
        cost_function_names{6} = 'volume_bounding_box'; % 
        cost_function_names{7} = 'weight_effort'; % 
        cost_function_names{8} = 'time_effort'; % 
        cost_function_names{9} = 'space_effort'; % 
        cost_function_names{10} = 'x_anchor_1'; % Distance to pick location
        cost_function_names{11} = 'rot_anchor_1'; % Distance to pick hand orientation
        cost_function_names{12} = 'x_anchor_2'; % Distance to placement location 
        cost_function_names{13} = 'rot_anchor_2'; % Distance to placement hand orientation
        cost_function_names{14} = 'com'; % 
        cost_function_names{15} = 'x_displace'; %
        cost_function_names{16} = 'dx_displace'; %
        cost_function_names{17} = 'dcom'; % 
     
        ccost_array{1} = ones(length(cost_function_names), 1)/length(cost_function_names);    
        
     case '22'
         
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'dddx';
        cost_function_names{4} = 'tau';
        cost_function_names{5} = 'dtau';
        cost_function_names{6} = 'ddtau'; %'';
        cost_function_names{7} = 'dq_tau'; %
        cost_function_names{8} = 'kinetic_en_sqrt'; % 
        cost_function_names{9} = 'quantity_motion_dx'; % 
        cost_function_names{10} = 'volume_bounding_box'; % 
        cost_function_names{11} = 'weight_effort'; % 
        cost_function_names{12} = 'time_effort'; % 
        cost_function_names{13} = 'space_effort'; % 
        cost_function_names{14} = 'flow_effort'; % 
     
        cost_function_names{15} = 'x_anchor_1'; % 
        cost_function_names{16} = 'x_anchor_2'; % 
        cost_function_names{17} = 'com'; % 
     
        cost_function_names{18} = 'x_displace'; %
        cost_function_names{19} = 'dx_displace'; %
        cost_function_names{20} = 'cartCurv'; % 
        cost_function_names{21} = 'dcom'; % 
        cost_function_names{22} = 'shapeDir'; % 
     
        ccost_array{1} = ones(length(cost_function_names), 1)/length(cost_function_names);
        
     case '26'
         
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'ddx';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'x_cartCurv';
        cost_function_names{6} = 'tau';
        cost_function_names{7} = 'dtau';
        cost_function_names{8} = 'ddtau'; %'';
        cost_function_names{9} = 'dq_tau'; %
        cost_function_names{10} = 'kinetic_en_sqrt'; % 
        cost_function_names{11} = 'quantity_motion_dx'; % 
        cost_function_names{12} = 'volume_bounding_box'; % 
        cost_function_names{13} = 'weight_effort'; % 
        cost_function_names{14} = 'time_effort'; % 
        cost_function_names{15} = 'space_effort'; % 
        cost_function_names{16} = 'flow_effort'; % 
        cost_function_names{17} = 'x_anchor_1'; % Distance to pick location
        cost_function_names{18} = 'rot_anchor_1'; % Distance to pick hand orientation
        cost_function_names{19} = 'x_anchor_2'; % Distance to placement location 
        cost_function_names{20} = 'rot_anchor_2'; % Distance to placement hand orientation
        cost_function_names{21} = 'com'; % 
     
        cost_function_names{22} = 'x_displace'; %
        cost_function_names{23} = 'dx_displace'; %
        cost_function_names{24} = 'cartCurv'; % 
        cost_function_names{25} = 'dcom'; % 
        cost_function_names{26} = 'shapeDir'; % 
     
        ccost_array{1} = ones(length(cost_function_names), 1)/length(cost_function_names);      
end

variableFactors.outputString = setting_outputstring;
variableFactors.win_length = 61; % 61 121 181 -- 101 201
variableFactors.win_shift = 20;
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

switch setting_dataset
    case 'sim'
        specStruct.dataset = 'Healthy1Ankle';
        specStruct.patient = [2];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 1;
        runMode = 'sim';
        
    case {'iit'}
        specStruct.dataset = 'IIT_2017';
        specStruct.patient = [1:2];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 1;
        runMode = 'win';
        
    case {'iit2d'}
        specStruct.dataset = 'IIT_2017_2D';
        specStruct.patient = [2];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 1;
        runMode = 'win';        
        
    case 'win2015' % non fatigue squats
        specStruct.dataset = 'Squats_URFI_2011';
        specStruct.patient = [1];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 0;
        runMode = 'win';
        
    case 'jump2D' % jumping to a target dataset
        specStruct.dataset = 'jump2D';
        specStruct.patient = [1];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 0;
        runMode = 'win';
        
    case 'expressive_ioc' % jumping to a target dataset
        specStruct.dataset = 'expressive_ioc';
        specStruct.patient = [1];
        specStruct.exerciseAcceptPrefixString = {''};
        specStruct.exerciseAcceptPrefixInstance = 0;
        runMode = 'win';        
end

variableFactors.datasetTag = specStruct.dataset;

variableFactors.normalizationMethod_feature = 'none';
variableFactors.normalizationMethod_cf = 'rset';
   
variableFactors.pivotMethod = setting_iocvalidation;
variableFactors.pivotSelection = 'all'; 
variableFactors.optThreshold = deg2rad(1e-6); % corresponds to deg2rad(1e-6)

variableFactors.simDistance = 100;
variableFactors.dynamicSampleRate = 0;

variableFactors.batchWindowBreakdownFactor = 3000; % (3000) break down the batch runs to increments of this factor, \pm 2 seconds on both sides
variableFactors.batchWindowBreakdownFactorOverlap = 600;

main_batch(specStruct, runMode, cost_function_names, setting_outputstring, variableFactors, setting_pathToRawData, setting_pathToOutput);
