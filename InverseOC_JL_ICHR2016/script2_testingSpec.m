function [specStruct, variableFactors] = script2_testingSpec(outputString, wantVariableOnly)

addpath(genpath(fullfile('..','Common')));

variableFactors.outputString = outputString;
outputStrSplit = strsplit(outputString, '_');

switch outputStrSplit{3}
    case '10'
        docsim_gen_settings_10knots; % 10knots 20knots
        variableFactors.knots = 10; % 4 7 10 12 // 10 20
        
    case '20'
        docsim_gen_settings_20knots; % 10knots 20knots
        variableFactors.knots = 20; % 4 7 10 12 // 10 20
        
    otherwise
        docsim_gen_settings_20knots;
        variableFactors.knots_sim = str2num(outputStrSplit{3});
        variableFactors.knots = str2num(outputStrSplit{3}); % 4 7 10 12 // 10 20
end

cconst_array = [1 1 1 1 1 1];

switch outputStrSplit{7}
    case '3a'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'ddx';
        cost_function_names{3} = 'tau';
        
        ccost_array{1} = [1 0 0];
        ccost_array{2} = [0 1 0];
        ccost_array{3} = [0 0 1];
        ccost_array{4} = [1 1 1];
        
%         cconst_array = [1 1 1 1 1 1];

    case '3b'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'ddx';
        cost_function_names{3} = 'tau';
        
        ccost_array{1} = [1 1 1]; 
        
    case '3c'
        cost_function_names{1} = 'dddq';
        cost_function_names{2} = 'ddx';
        cost_function_names{3} = 'tau';
        
        ccost_array{1} = [1 0 0];  
        
    case '3e'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'ddx';
        cost_function_names{3} = 'tau';
        
        ccost_array{1} = [0 1 0];
        
    case '3f'
        cost_function_names{1} = 'dddq';
        cost_function_names{2} = 'dddx';
        cost_function_names{3} = 'tau';
        
        ccost_array{1} = [1 0 0];
        ccost_array{2} = [0 1 0];
        ccost_array{3} = [0 0 1];   
                ccost_array{4} = [1 1 1];
        
    case '3g'
        cost_function_names{1} = 'tau';
        cost_function_names{2} = 'dtau';
        cost_function_names{3} = 'ddtau';
        
        ccost_array{1} = [1 0 0];
        ccost_array{2} = [0 1 0];
        ccost_array{3} = [0 0 1];     
                ccost_array{4} = [1 1 1];
    case '3h'
        cost_function_names{1} = 'dq_tau';
        cost_function_names{2} = 'kinetic_en_sqrt';
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
        
    case {'9', '9a'}
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'ddx';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'tau';
        cost_function_names{6} = 'dtau';
        cost_function_names{7} = 'ddtau'; %'';
        cost_function_names{8} = 'dq_tau'; %
        cost_function_names{9} = 'kinetic_en_sqrt'; % 
     
        ccost_array{1} = [1 0 0 0 0 0 0 0 0];  
        
    case '9b'
        cost_function_names{1} = 'ddx';
        cost_function_names{2} = 'dq_tau'; %
        cost_function_names{3} = 'kinetic_en_sqrt'; % 
        cost_function_names{4} = 'ddq';
        cost_function_names{5} = 'dddq';        
        cost_function_names{6} = 'dddx';
        cost_function_names{7} = 'tau';
        cost_function_names{8} = 'dtau';
        cost_function_names{9} = 'ddtau'; %'';
       
     
        ccost_array{1} = [1 0 0 0 0 0 0 0 0];  
        
    case '9c'
        cost_function_names{1} = 'tau';
        cost_function_names{2} = 'dtau';
        cost_function_names{3} = 'ddtau'; %'';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'ddq';
        cost_function_names{6} = 'dddq';
        cost_function_names{7} = 'ddx';
        cost_function_names{8} = 'dq_tau'; %
        cost_function_names{9} = 'kinetic_en_sqrt'; %
        
        ccost_array{1} = [1 0 0 0 0 0 0 0 0];    
        
    case '9I'
        cost_function_names{1} = 'tau';
        cost_function_names{2} = 'dtau';
        cost_function_names{3} = 'ddtau'; %'';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'ddq';
        cost_function_names{6} = 'dddq';
        cost_function_names{7} = 'ddx';
        cost_function_names{8} = 'dq_tau_sq'; %
        cost_function_names{9} = 'kinetic_en_sqrt_sq'; %
        
        ccost_array{1} = [1 0 0 0 0 0 0 0 0];          
        
    case '11'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'ddx';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'tau';
        cost_function_names{6} = 'dtau';
        cost_function_names{7} = 'ddtau'; %'';
        cost_function_names{8} = 'dq_tau'; %  
        cost_function_names{9} = 'kinetic_en_sqrt'; %
        cost_function_names{10} = 'cop'; %
        cost_function_names{11} = 'dcop'; %
        
    case '12'
        cost_function_names{1} = 'ddq';
        cost_function_names{2} = 'dddq';
        cost_function_names{3} = 'ddx';
        cost_function_names{4} = 'dddx';
        cost_function_names{5} = 'tau';
        cost_function_names{6} = 'dtau';
        cost_function_names{7} = 'ddtau'; %'';
        cost_function_names{8} = 'dq_tau'; %  
        cost_function_names{9} = 'kinetic_en_sqrt'; %
        cost_function_names{10} = 'potential_en'; %
        cost_function_names{11} = 'cop'; %
        cost_function_names{12} = 'dcop'; %
        
        ccost_array{1} = [1 0 0 0 0  0 0 0 0 0 0 0];           
end

variableFactors.win_length = str2num(outputStrSplit{4}); % 61 121 181 -- 101 201
variableFactors.win_shift = str2num(outputStrSplit{5});
if variableFactors.win_shift == 0
    variableFactors.win_shift = floor(variableFactors.win_length/2);
end
variableFactors.ioc_cf_q_constraint = 'startmid1endonknot'; %  (fixed, variable, edge, startmidend, startmidmidend) where to set the consts? 3 points (fixed), global matching constraint (variable).

if strcmpi(variableFactors.ioc_cf_q_constraint, 'sm1ek')
    variableFactors.ioc_cf_q_constraint = 'startmid1endonknot';
end

variableFactors.generateSimOnly = 0;
variableFactors.spline_length = variableFactors.win_length;

variableFactors.correlation_threshold = str2num(outputStrSplit{9})/10000; % [0:1]
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

variableFactors.segment_only_windows = outputStrSplit{8}; % (none, constraint, seg, mergeseg) if yes, then win_length/spline_length/win_shift is ignored and adjusted according to some manual segment
variableFactors.windowMode = 'wide'; % wide narrow

specStruct.ccost_array = ccost_array;
specStruct.cconst_array = cconst_array;


outputStrSplit_dataset = strsplit(outputStrSplit{1}, 'q');

switch outputStrSplit_dataset{1}
    case 'sim'
        specStruct.dataset = 'Squats_TUAT_2015';
        specStruct.patient = [1];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 1;
        runMode = 'sim';
        
    case 'iit'
        specStruct.dataset = 'IIT_2017';
        specStruct.patient = [2];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 1;
        runMode = 'win';
        
    case 'win2015' % non fatigue squats
        specStruct.dataset = 'Squats_TUAT_2015';
        specStruct.patient = 1:8;
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 0;
        runMode = 'win';
        
    case 'win2011' % fatigue squats
        specStruct.dataset = 'Squats_TUAT_2011';
        specStruct.patient = 2:7;
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 0;
        runMode = 'win';
               
    case 'winh1ankle' % healthy1
        specStruct.dataset = 'Healthy1Ankle';
        specStruct.patient = 1:20;
%         specStruct.exerciseAcceptPrefixString = {'HFEO_SUP_SLO', 'KEFO_SIT_SLO', 'KHEF_SUP_SLO', 'SQUA_STD_SLO', 'STSO_SIT_SLO'};
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD_SLO', 'STSO_SIT_SLO'};
        specStruct.exerciseAcceptPrefixInstance = [0 0 0 0 0];
        runMode = 'win';     
                
    case 'winh1ankle2' % healthy1
        specStruct.dataset = 'Healthy1Ankle';
        specStruct.patient = 1;
%         specStruct.exerciseAcceptPrefixString = {'HFEO_SUP_SLO', 'KEFO_SIT_SLO', 'KHEF_SUP_SLO', 'SQUA_STD_SLO', 'STSO_SIT_SLO'};
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD_SLO'};
        specStruct.exerciseAcceptPrefixInstance = [1];
        runMode = 'win';  
         
    case 'winh1hip' % healthy1
        specStruct.dataset = 'Healthy1Hip';
        specStruct.patient = 1:20;
        specStruct.exerciseAcceptPrefixString = {'HFEO_SUP_SLO', 'KEFO_SIT_SLO', 'KHEF_SUP_SLO'};
        specStruct.exerciseAcceptPrefixInstance = [0 0 0 0 0];
        runMode = 'win';      
        
    case 'winh1hip1' % healthy1
        specStruct.dataset = 'Healthy1Hip';
        specStruct.patient = 3;
        specStruct.exerciseAcceptPrefixString = {'HFEO_SUP_SLO', 'KEFO_SIT_SLO', 'KHEF_SUP_SLO'};
        specStruct.exerciseAcceptPrefixInstance = [0 0 0 0 0];
        runMode = 'win';
        
    case 'winh1hip2' % healthy1
        specStruct.dataset = 'Healthy1Hip';
        specStruct.patient = 3;
        specStruct.exerciseAcceptPrefixString = {'HFEO_SUP_SLO', 'KEFO_SIT_SLO', 'KHEF_SUP_SLO'};
        specStruct.exerciseAcceptPrefixInstance = [0 0 0 0 0];
        runMode = 'win';
        
    case 'winh1hip3' % healthy1
        specStruct.dataset = 'Healthy1Hip';
        specStruct.patient = 3;
        specStruct.exerciseAcceptPrefixString = {'HFEO_SUP_SLO', 'KEFO_SIT_SLO', 'KHEF_SUP_SLO'};
        specStruct.exerciseAcceptPrefixInstance = [0 0 0 0 0];
        runMode = 'win';      
end

variableFactors.datasetTag = outputStrSplit_dataset{1};

switch outputStrSplit_dataset{1}(1:3)
    case 'win'
        if length(outputStrSplit_dataset) == 1
%             specStruct.patient = [1];
        else
            newPxArray = [];
%             if length(outputStrSplit_dataset{2}) == 2
%                 for ind_px = 1:length(outputStrSplit_dataset{2})
%                     newPxArray = [newPxArray str2num(outputStrSplit_dataset{2}(ind_px))];
%                 end
%             else
                for ind_px = 2:length(outputStrSplit_dataset)
                    newPxArray = [newPxArray str2num(outputStrSplit_dataset{ind_px})];
                end
%             end

            specStruct.patient = newPxArray;
        end
end

if length(outputStrSplit{2}) == 1
    outputStrSplit{2} = [outputStrSplit{2} 'r'];
end

switch specStruct.dataset
    case 'Squats_TUAT_2011'
        switch outputStrSplit{2}(1)
            case '1'
                variableFactors.dofsToUse = 3; % 3, 2:4
                
            case '2' % ankle, knee
                variableFactors.dofsToUse = 2:3; % 3, 2:4
                
            case '3' % ankle, knee, hip
                variableFactors.dofsToUse = 2:4; % 3, 2:4
                
            case '4'
                variableFactors.dofsToUse = 1:4;
                
            case '5' % ankle, knee, hip, torso (pelvis is not a joint in 2011)
                variableFactors.dofsToUse = 2:6; % 3, 2:4
                
            case '6' % ankle, knee, hip, torso, shoulder (pelvis is not a joint)
                variableFactors.dofsToUse = 2:7; % 3, 2:4
        end
        
    case 'Squats_TUAT_2015'
        switch outputStrSplit{2}(1)
            case '1'
                variableFactors.dofsToUse = 3; % 3, 2:4
                
            case '2' % ankle, knee
                variableFactors.dofsToUse = 2:3; % 3, 2:4
                
            case '3' % ankle, knee, hip
                variableFactors.dofsToUse = 2:4; % 3, 2:4
                
            case '4'
                variableFactors.dofsToUse = 1:4;
                
            case '5' % ankle, knee, hip, pelvis, torso 
                variableFactors.dofsToUse = 2:6; % 3, 2:4
                
            case '6' % ankle, knee, hip, pelvis, torso, shoulder 
                variableFactors.dofsToUse = 2:7; % 3, 2:4
        end
        
    case {'healthy1hip', 'healthy1ankle'}
        switch outputStrSplit{2}(1)
            case '3' % ankle, knee, hip
                variableFactors.dofsToUse = 2:4; % 3, 2:4
                
            case '4'  
                variableFactors.dofsToUse = 1:4; 
        end
        
    case 'IIT_2017'
        variableFactors.dofsToUse = 7:45;
end

variableFactors.normalizationMethod_feature = outputStrSplit{10};
variableFactors.normalizationMethod_cf = outputStrSplit{11};
    
switch outputStrSplit{2}(2)

    case 'm'
         variableFactors.pivotMethod = 'resnorm_rmse1'; % resnorm rmse
         variableFactors.dataLoadLength = 0;
         variableFactors.pivotSelection = 'all'; 
         variableFactors.optThreshold = deg2rad(1e-6); % corresponds to deg2rad(1e-6)

     case 'r'
         variableFactors.pivotMethod = 'resnorm'; % resnorm rmse
         variableFactors.dataLoadLength = 0;
         variableFactors.pivotSelection = 'all'; 
         variableFactors.optThreshold = deg2rad(1e-6); % corresponds to deg2rad(1e-6)
end

variableFactors.simDistance = str2num(outputStrSplit{12});
variableFactors.dynamicSampleRate = str2num(outputStrSplit{13});

variableFactors.batchWindowBreakdownFactor = 3000; % (3000) break down the batch runs to increments of this factor, \pm 2 seconds on both sides
variableFactors.batchWindowBreakdownFactorOverlap = 600;

% specStruct.exerciseAcceptPrefixInstance = [1 1 1 1 1];

if ~exist('wantVariableOnly', 'var')
    mainDataFolder = 'D:\aslab\data';
    mainOutputFolder = 'D:\results\IOC\IOC51_modelRL_IIT2';
    main_batch(specStruct, runMode, cost_function_names, outputString, variableFactors, mainDataFolder, mainOutputFolder);
end