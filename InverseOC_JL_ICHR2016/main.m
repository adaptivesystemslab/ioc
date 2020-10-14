function main(filesToLoad, manSegLoadPath, outputPath, currInstName, runSettings, currFilestack, run_mode, cost_function_names)

% to modify the cost function, update 'cost_function_names', 'ccost_array', 
% and 'J_array' in calc_direct_cost, and it should automatically proprogate
% however, if simulation is being used, the output plots will not
% automatically update and needs to be changed by hand

safetySave = 1e10; % after this much interations, dump the memory to file instead of only waiting until the end
saveFolder = 'saved';
checkMkdir(saveFolder); 

runFirstDOCOnly = 0;
shortWinThreshold = 5;

if ~exist('cost_function_names', 'var')
    cost_function_names{1} = 'dddx';
    cost_function_names{2} = 'dddq';
    cost_function_names{3} = 'ddq';
    cost_function_names{4} = 'tau';
    cost_function_names{5} = 'dtau';
    cost_function_names{6} = 'geodesic';
    cost_function_names{7} = 'energy';
    cost_function_names{8} = 'effort';
    % cost_function_names{9} = 'dcop';
end

feature_load_names{1} = 't';
feature_load_names{2} = 'q';
feature_load_names{3} = 'dq';
feature_load_names{4} = 'ddq';
feature_load_names{5} = 'dddq';
feature_load_names{6} = 't_cf_q_const';
feature_load_names{7} = 't_cf_dq_const';
feature_load_names{8} = 't_cf_ddq_const';
feature_load_names{9} = 't_knot';
feature_load_names{10} = 't_sp_dq_const';
feature_load_names{11} = 'y_cf_q_const';
feature_load_names{12} = 'y_cf_dq_const';
feature_load_names{13} = 'y_cf_ddq_const';
feature_load_names{14} = 'y_knot';
feature_load_names{15} = 'y_sp_dq_const';

% rmse_fct = @(a,b,c) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2)))*(100/c); % a, b: the two data arrays, c: pre-spline length
rmse_fct = @(a,b,c) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2))); % a, b: the two data arrays, c: pre-spline length
% nrmse_fct = @(a,b) 100*(sqrt(sum((a-b).^2)/length(a)))/sqrt(sum((a).^2)/length(a));
% CC_fct = @(a,b) corr2(a,b);

diary(fullfile(outputPath, [currInstName '_diary.txt'])); % write console to file
outputPath_intermedFig = fullfile(outputPath, 'intermed_fig');
outputPath_intermedCsv = fullfile(outputPath, 'intermed_csv');
outputPath_intermedMat = fullfile(outputPath, 'intermed_mat');
checkMkdir(outputPath_intermedFig);
checkMkdir(outputPath_intermedCsv);
checkMkdir(outputPath_intermedMat);

intermedMatPath = {};

% load the initial parameters
runSettings_docsim = runSettings; % copy out the sim settings before
switch runSettings.variableFactors.half_phase_construction
    case 'yes'
        runSettings_docsim.variableFactors.knots = ceil(runSettings.variableFactors.knots/2);
        
%     case 'double'
%         runSettings_docsim.variableFactors.knots = ceil(runSettings.variableFactors.knots*2);
        
    case 'no'
        
end

cconst = currFilestack.cconst_array;  

[param, traj_load, segmentInfo] = setup_main(filesToLoad, manSegLoadPath, currFilestack, run_mode, runSettings, cost_function_names, feature_load_names);
% param.cost_functions_ioc = 1:length(cost_function_names); % this initalizes the cost functions to be used. the IOC may remove certain cost functions due to colinearity, and this variables tracks that

runSettings.variableFactors
    
% create a separate param file to use here, and apply any local
% modifications needed
runSettings_docsim.variableFactors.doc_sim_win_length = runSettings.variableFactors.doc_sim_win_length_sim;
runSettings_docsim.variableFactors.spline_length = runSettings.variableFactors.spline_length_sim;
runSettings_docsim.variableFactors.win_length = runSettings.variableFactors.win_length_sim;
[param_docsim, ~, ~] = setup_main(filesToLoad, manSegLoadPath, currFilestack, run_mode, runSettings_docsim, cost_function_names, feature_load_names);
[param_docsim, ccost_array, const_x_array, const_y_array] = halfPhaseConstruction(currFilestack, param_docsim);
param_docsim.coeff_cf = param.coeff_cf; % the simulation normalization can vary based on the length of window being passed in (ie param.win_length)

switch run_mode
    case 'sim'
        % set the cost functions for simulation testing: ddq, ddx, tau
        
        ccost_runningStrSum = 0;
        for ind_ccost = 1:length(currFilestack.ccost_array)
            ccost_runningStrSum = ccost_runningStrSum + sum(currFilestack.ccost_array{ind_ccost});
        end
        
%         create the savefile string name
        simStr = ['L' num2str(param_docsim.doc_sim_win_length) ...
            '_CF' arrayToStr(length(cost_function_names)) ...
            '_CFN' cost_function_names{1} '_' cost_function_names{2} '_' cost_function_names{end} ...
            '_CFC' num2str(ccost_array{1}(1)) '_' num2str(ccost_array{1}(2)) '_' num2str(ccost_array{1}(3)) ...
            '_' param.normalizationMethod_feature '_' param.normalizationMethod_cf '_' param.outputString];
   
%         simStr = 'L2.01_CF3_CFNddq_ddsx_tau_CFC1_0_0_none_rset_sim_3r_10_61_20_startmid1endonknot_3a_none_10000_none_rset_100_00_20170215092359';
%         simStr = 'L2.01_CF3_CFNddq_ddx_tau_CFC1_0_0_none_rset_sim_3r_10_61_20_startmid1endonknot_3a_none_10000_none_rset_100_10_20170201020725';
%         simStr = 'L2.01_CF3_CFNdddq_ddx_tau_CFC1_0_0_none_rset_20170131192953sim_3r_20_61_20_startmid1endonknot_3c_none_10000_none_rset_100_0';
%         simStr = 'L2.01_CF3_CFNdddq_ddx_tau_CFC1_0_0_none_rset_20170131202800sim_3r_20_61_20_startmid1endonknot_3c_none_10000_none_rset_100_5';
                
        docsim_path = fullfile(saveFolder, [simStr '.mat']);
            
%         if strcmpi(computer, 'PCWIN')
%             % only run this if it's not on sharcnet
% %             param_docsim_hash = DataHash(param_docsim);
% %             docsim_path = fullfile(saveFolder, [param_docsim_hash '.mat']);
% docsim_path = fullfile(saveFolder, [param_docsim_hash '.mat']);
%             simStr
%         else
%             docsim_path = [];
%         end
        
        if exist(docsim_path, 'file') 
            fprintf('File profile match database, loading %s\n', docsim_path);
            load(docsim_path) % load file if exist
            
            for ind_la = 1:length(set_docsim_check_struct)
                set_docsim_check_struct{ind_la}
            end
        else
            % otherwise, generate the trajectory
        fprintf('File profile not found, generating %s\n', docsim_path);
            
        % construct the sequence of optimal trajectory based on the ccost_array
        windowCount = length(ccost_array);
        t_full = [];
        t_cf_q_const = [];
        segmentInfo.timeStart = [];
        segmentInfo.timeEnd = [];
        
        for ind_windowCount = 1:windowCount
            ccost_curr = ccost_array{ind_windowCount};
            param_docsim.const_x = const_x_array{ind_windowCount};
            param_docsim.const_y = const_y_array{ind_windowCount};
            
            fprintf('(%u/%u) Solving the direct problem...\n', ind_windowCount, windowCount);
            
            [param_docsim, set_docsim_check_struct{ind_windowCount}] = set_docsim_constraints(param_docsim);
            feature_win = [];
            
            [feature_optsim, J_array_optsim, J_optsim_contrib, output_optsim(ind_windowCount)] = main_direct(ccost_curr, cconst, feature_win, param_docsim, param_docsim.doc_sim);
            J_optsim_contrib_array{ind_windowCount} = J_optsim_contrib;
            J_optsim_array_array{ind_windowCount} = J_array_optsim';

            % upsample the function and include blanks (if needed)
            newSampleRate = param_docsim.dt_spline;
            lengthBlank = 0; 
            
            t_opt_upsample = param_docsim.ti_full:newSampleRate:param_docsim.tf_full;
            t_opt_blanks = (1:lengthBlank)*newSampleRate + t_opt_upsample(end);
%             q_opt_upsamp = interp1(param.t, q_optsim', t_opt_upsample, 'spline')';            
%             q_opt_blanks = repmat(q_opt_upsamp(:, end), [1, lengthBlank]);
%             q_full_cell{ind_windowCount_gen} = [q_optsim q_opt_blanks]; % save each generated traj
%             t_opt_use = [t_opt_upsample t_opt_blanks]; 
            
            q_full_cell{ind_windowCount} = feature_optsim.q;
            dq_full_cell{ind_windowCount} = feature_optsim.dq;
            ddq_full_cell{ind_windowCount} = feature_optsim.ddq;
            dddq_full_cell{ind_windowCount} = feature_optsim.dddq;

            t_knot_ind{ind_windowCount} = feature_optsim.xknot_ind; 
            
            % set up the time array
            if ind_windowCount == 1
                t_offset = 0;
                t_opt_upsample_len = length(t_opt_upsample);
                t_opt_blanks_len = length(t_opt_blanks);
            elseif 0
                % if want no overlap with previous window
                t_offset = t_full_cell{end}(end)+param_docsim.dt_spline;
                t_opt_upsample_len = length(t_opt_upsample);
                t_opt_blanks_len = length(t_opt_blanks);
            else
                % if want overlap with previous window
                t_offset = t_full_cell{end}(end); % the last -1 induces an overlap of 1 step to the previous window
                t_opt_upsample_len = length(t_opt_upsample) - 1;
                t_opt_blanks_len = length(t_opt_blanks) - 1;
            end
            
            t_full_cell{ind_windowCount} =    t_offset + param_docsim.t_spline;
            t_knot{ind_windowCount} =         t_offset + feature_optsim.xknot_t;
            t_cf_q_const{ind_windowCount} =   t_offset + param_docsim.t_spline(param_docsim.const_x);
            t_cf_dq_const{ind_windowCount} =  t_offset + param_docsim.t_spline(param_docsim.const_dx);
            t_cf_ddq_const{ind_windowCount} = t_offset + param_docsim.t_spline(param_docsim.const_ddx);
            y_cf_q_const{ind_windowCount} =   param_docsim.const_y;
            y_cf_dq_const{ind_windowCount} =  param_docsim.const_dy;
            y_cf_ddq_const{ind_windowCount} = param_docsim.const_ddy;
            t_sp_dq_const{ind_windowCount} =  t_offset + param_docsim.t_spline(param_docsim.spline_endcond_dx);
            
            % now create an array for the ground truth c_cost and
            % contribution ratio for plotting later on, and assessing error
            ccost_opt_upsample = repmat(ccost_curr',                  [1 t_opt_upsample_len]);
            ccost_opt_blanks =   repmat(ones(length(ccost_curr), 1),  [1 t_opt_blanks_len]);
            ccost_opt_cell{ind_windowCount} = [ccost_opt_upsample ccost_opt_blanks];
            
            J_opt_upsample = repmat(J_optsim_contrib',                [1 t_opt_upsample_len]);
            J_opt_blanks =   repmat(zeros(length(ccost_curr), 1),     [1 t_opt_blanks_len]);
            J_opt_cell{ind_windowCount} = [J_opt_upsample J_opt_blanks];
            
            output_optsim_mapping_upsample = repmat(ind_windowCount,  [1 t_opt_upsample_len]);
            output_optsim_mapping_blanks =   repmat(ind_windowCount,  [1 t_opt_blanks_len]);
            output_optsim_mapping{ind_windowCount} = [output_optsim_mapping_upsample output_optsim_mapping_blanks];
            
            % track segment windows for 'param.segment_only_widows' mode
%             segmentInfo.timeStart(ind_windowCount) = t_full(end-length(param_docsim.t_spline)+1);
%             segmentInfo.timeEnd(ind_windowCount) = t_full(end);
            segmentInfo.timeStart(ind_windowCount) = t_full_cell{ind_windowCount}(1);
            segmentInfo.timeEnd(ind_windowCount) = t_full_cell{ind_windowCount}(end);
        end
        
        % combine all the generated trajectory into a single var
        t_full = horzcat(t_full_cell{:});
        q_full = horzcat(q_full_cell{:});
        dq_full = horzcat(dq_full_cell{:});
        ddq_full = horzcat(ddq_full_cell{:});
        dddq_full = horzcat(dddq_full_cell{:});
%         dq_full = calcDeriv(q_full, param_docsim.dt_spline);
%         ddq_full = calcDeriv(dq_full, param_docsim.dt_spline);
         
        t_knot_full = horzcat(t_knot{:});
        t_q_const_full = horzcat(t_cf_q_const{:});
        t_dq_const_full = horzcat(t_cf_dq_const{:});
        t_ddq_const_full = horzcat(t_cf_ddq_const{:});
        t_endcond_full = horzcat(t_sp_dq_const{:});
        
        y_q_const_full = horzcat(y_cf_q_const{:});
        y_dq_const_full = horzcat(y_cf_dq_const{:}); 
        y_ddq_const_full = horzcat(y_cf_ddq_const{:}); 
        
        % check for unique entries
        [t_full, t_full_IA] = unique(t_full);
        [t_knot_full, t_knot_full_IA] = unique(t_knot_full);
        [t_q_const_full, t_q_const_full_IA] = unique(t_q_const_full);
        [t_dq_const_full, t_dq_const_full_IA] = unique(t_dq_const_full);
        [t_ddq_const_full, t_ddq_const_full_IA] = unique(t_ddq_const_full);
        [t_endcond_full, t_endcond_full_IA] = unique(t_endcond_full);
        
        q_full = q_full(:, t_full_IA);
        dq_full = dq_full(:, t_full_IA);
        ddq_full = ddq_full(:, t_full_IA);
        dddq_full = dddq_full(:, t_full_IA);
        y_q_const_full = y_q_const_full(:, t_q_const_full_IA);
        y_dq_const_full = y_dq_const_full(:, t_dq_const_full_IA);
        y_ddq_const_full = y_ddq_const_full(:, t_ddq_const_full_IA);
        
        traj_load_lengthAdj.t = t_full;
        traj_load_lengthAdj.q = q_full;
        traj_load_lengthAdj.dq = dq_full;
        traj_load_lengthAdj.ddq = ddq_full;
        traj_load_lengthAdj.dddq = dddq_full;
        traj_load_lengthAdj.t_knot = t_knot_full;
        traj_load_lengthAdj.t_cf_q_const = t_q_const_full;
        traj_load_lengthAdj.t_cf_dq_const = t_dq_const_full;
        traj_load_lengthAdj.t_cf_ddq_const = t_ddq_const_full;
        traj_load_lengthAdj.t_sp_dq_const = t_endcond_full;
        traj_load_lengthAdj.y_knot = [];
        traj_load_lengthAdj.y_cf_q_const = y_q_const_full;
        traj_load_lengthAdj.y_cf_dq_const = y_dq_const_full;
        traj_load_lengthAdj.y_cf_ddq_const = y_ddq_const_full;
        traj_load_lengthAdj.y_sp_dq_const = [];  
        
        ccost_full = horzcat(ccost_opt_cell{:})'; 
        J_full = horzcat(J_opt_cell{:})';
        output_optsim_mapping_full = horzcat(output_optsim_mapping{:})'; 
        
        if ~isempty(docsim_path)
            fprintf('Saving template at %s \n', docsim_path); 
            save(docsim_path, 'traj_load_lengthAdj', 'ccost_full', 'J_full', 'output_optsim_mapping_full', 'set_docsim_check_struct', ...
                'output_optsim', 'segmentInfo', 'traj_load_lengthAdj', 't_full_cell', 'q_full_cell', 'dq_full_cell', 't_knot_ind', ...
                't_q_const_full', 'y_q_const_full', 't_dq_const_full', 'y_dq_const_full', 'J_optsim_contrib_array', 'J_optsim_array_array');
        end
        end
        
    case 'win'        
        % load the full trajectory
%         traj_load = load(trajToLoadPath);
        
        lengthToKeep = param.dataLoadLength; 
        if isempty(lengthToKeep)
            lengthToKeep = 1:length(param.t_full);
        end 
        
        % generate the t_knot based on the simulation settings
        [param_docsim, skip_count] = update_intermed_ind(param_docsim, param_docsim.doc_sim.knot_count);
        length_traj = length(param.t_full);
        intermed_ind = unique([1:skip_count:length_traj-1 length_traj]); % this ensures there is a knot point at the end, and on constraints
        
        xknot_t = param.t_full(intermed_ind);
        
        traj_load_lengthAdj.t = param.t_full(lengthToKeep);
        traj_load_lengthAdj.q = traj_load.q(:, lengthToKeep);
        traj_load_lengthAdj.dq = traj_load.dq(:, lengthToKeep);
        traj_load_lengthAdj.ddq = traj_load.ddq(:, lengthToKeep);
        traj_load_lengthAdj.dddq = traj_load.dddq(:, lengthToKeep);
        traj_load_lengthAdj.t_knot = xknot_t;
        
        % remove segment windows that don't fall into the spec'ed area
%         segmentInfo = removeOverflowSegments(traj_load_lengthAdj.t, segmentInfo);

        ccost_full = ones(length(cost_function_names), length(traj_load_lengthAdj.t))';
        J_full = zeros(length(cost_function_names), length(traj_load_lengthAdj.t))';
        
        output_optsim_mapping_full = ones(length(cost_function_names), length(traj_load_lengthAdj.t))';
end

if isfield(runSettings.variableFactors, 'generateSimOnly') && runSettings.variableFactors.generateSimOnly == 1
    % this function is triggered to only run and save the test data
    return
end

feature_full = feature_collate(traj_load_lengthAdj, feature_load_names);

if 0 
    figure; plot(feature_full.t, feature_full.q');
    hold on
    plot(feature_full.t_cf_q_const, zeros(size(feature_full.t_cf_q_const)), 'x');
end

% convert to constraint point tracking
% param.segment_only_windows = 'seg';

% break it down into windows and figure out how much data to feed in at a
% time for the windowing approach
switch param.segment_only_windows
    case 'none'
        % create windows that are regular intervals and shifts
        ind = 0;
        for ind_windowCount = 1:param.win_shift:length(feature_full.t)
            ind = ind+1;
            indToUse_window(ind, 1) = ind_windowCount;
            indToUse_window(ind, 2) = ind_windowCount + param.win_length - 1;
            
            if indToUse_window(ind, 2) > length(feature_full.t)
                %         indToUse_window(ind, 2) = length(t_full);
                indToUse_window = indToUse_window(1:ind-1, :); % remove that latest entry
                break % we're at the end of the line
            end
        end
        
    case 'fixed_constraints'
           % create windows that are regular intervals and shifts
     
           indToUse_window = [1	61
               21	81
               41	101
               61	121
               101	161
               121	181
               141	201
               161	221
               201	261
               221	281
               241	301
               261	321
               301	361
               321	381
               341	401];

    case 'constraint'
        % create windows that window by segment definition only
        ind_counter = 0;
        for ind_windowCount = 1:length(feature_full.t_cf_q_const)-1
            [startVal, startInd] = findClosestValue(feature_full.t_cf_q_const(ind_windowCount), feature_full.t);
            [endVal, endInd] = findClosestValue(feature_full.t_cf_q_const(ind_windowCount+1), feature_full.t);
            
            if endInd - startInd < shortWinThreshold
                % if two constraints are too close to each other
                continue
            end
            
            ind_counter = ind_counter + 1;
            indToUse_window(ind_counter, 1) = startInd;
            indToUse_window(ind_counter, 2) = endInd;
        end
        
    case 'seg'
        % create windows that window by segment definition only
        for ind_windowCount = 1:length(segmentInfo.timeStart)
            [startVal, startInd] = findClosestValue(segmentInfo.timeStart(ind_windowCount), feature_full.t);
            [endVal, endInd] = findClosestValue(segmentInfo.timeEnd(ind_windowCount), feature_full.t);

            indToUse_window(ind_windowCount, 1) = startInd;
            indToUse_window(ind_windowCount, 2) = endInd;
        end
        
   case 'mergeseg2'
        ind = 0;
        for ind_windowCount = 1:2:length(segmentInfo.timeStart)
            [startVal, startInd] = findClosestValue(segmentInfo.timeStart(ind_windowCount), feature_full.t);
            [endVal, endInd] = findClosestValue(segmentInfo.timeEnd(ind_windowCount+1), feature_full.t);

            ind = ind + 1;
            indToUse_window(ind, 1) = startInd;
            indToUse_window(ind, 2) = endInd;
        end        
        
   case 'mergeseg5'
        ind = 0;
        for ind_windowCount = 1:5:length(segmentInfo.timeStart)
            [startVal, startInd] = findClosestValue(segmentInfo.timeStart(ind_windowCount), feature_full.t);
            [endVal, endInd] = findClosestValue(segmentInfo.timeEnd(ind_windowCount+1), feature_full.t);

            ind = ind + 1;
            indToUse_window(ind, 1) = startInd;
            indToUse_window(ind, 2) = endInd;
        end   
        
    case 'knot_all'
        ind_counter = 0;
        for ind_windowCount = 1:length(feature_full.t_knot)-1
            [startVal, startInd] = findClosestValue(feature_full.t_knot(ind_windowCount), feature_full.t);
            
            for ind_windowCount2 = ind_windowCount+1:length(feature_full.t_knot)
                [endVal, endInd] = findClosestValue(feature_full.t_knot(ind_windowCount2), feature_full.t);
            
            if endInd - startInd < shortWinThreshold
                % if two constraints are too close to each other
                continue
            end
            
            ind_counter = ind_counter + 1;
            indToUse_window(ind_counter, 1) = startInd;
            indToUse_window(ind_counter, 2) = endInd;  
            end

        end
end

windowCount = size(indToUse_window, 1);

tic % start timing the IOC recovery process

% blank variables, for plotting later on
t_rmse = []; 
ccostRMSE_plot = [];
ccostAll_plot = [];
ratioRMSE_plot = [];
ratioAll_plot = [];
rmseMin_plot = [];
condMin_plot = [];
resnormMin_plot = [];
q_jh_absmax_plot = [];
rmseAll_plot = [];
condAll_plot = [];
resnormAll_plot = [];
dofsUsedInIOC_plot = [];
minRmseInd_plot = [];
q_jh_range_plot = [];

resnormAll_lsqlin_const_plot =   [];  
resnormAll_lsqlin_unconst_plot = []; 
resnormAll_lsqlin_const_RMSE_plot = [];
resnormAll_lsqlin_unconst_RMSE_plot = [];
c_ioc_constrain_plot = []; 
c_ioc_unconstrain_plot = []; 

param.feature_full = feature_full;

for ind_windowCount = 1:windowCount
    fprintf('(%u/%u) Loading window...\n', ind_windowCount, windowCount);
%     param = setup_main(trajToLoadPath, dynToLoadPath); % reset param in case something is remaining

    switch param.segment_only_windows
        case 'none'
            
        otherwise
            % update the length of the window
            new_win_len = indToUse_window(ind_windowCount, 2) - indToUse_window(ind_windowCount, 1) + 1;
            param = update_length(param, new_win_len);
    end

    % pull out the window of data
    feature_win = feature_windowed(feature_full, indToUse_window(ind_windowCount, 1), indToUse_window(ind_windowCount, 2), param);
    feature_win_save{ind_windowCount} = feature_win;
    
    docWinIndToUse = indToUse_window(ind_windowCount, 1) + 1; % this ensures that it's pulling the orig DOC info from the correct box
    ccost_win = ccost_full(docWinIndToUse, :);  
    
    if length(ccost_win) ~= length(cost_function_names)
        ccost_win = ones(size(cost_function_names));
    end
    
%     if min(feature_win.t) < 10
%         continue
%     end
    
%     if 0
        ccost_win_init = zeros(size(ccost_win));
%     elseif ind_windowCount == 1
%         use the global initial
%         ccost_win_init = ccost_win;
%     else
%         % use the weights from the previous window
%         ccost_win_init = ccostRMSE_select;
%     end

    % update the constraints based on windowed data
    [param, set_ioc_check_struct(ind_windowCount)] = set_ioc_constraints(feature_win, feature_full, indToUse_window(ind_windowCount, 1), param, 'ioc');   
    
%     param = update_intermed_ind(param, param.doc_pivot.knot_count); % if we're moving the constraint point, we also need to move the pivot
            
    if ind_windowCount == 1
        param.cost_functions_ioc = [];
    end

    % call the inverse problem
    saveprefix = fullfile(saveFolder, datestr(now, 'yyyymmddHHMMSS'));
    fprintf('Solving the inverse problem...\n');
    [ccost_inverse_pivot, J_inverse, output_inverse{ind_windowCount}, param.cost_functions_ioc, ...
         J_costfull_inverse{ind_windowCount}, J_costconst_check{ind_windowCount}] = ...
        main_inverse2(feature_win, ccost_win, cconst, param, param.ioc, ccost_win_init, [], saveprefix);
    
    
    
    switch param.pivotMethod
        case 'rmse'
            % call the direct problem
        direct_check_flags = {};
        [param, set_doc_rmse_check_struct(ind_windowCount)] = set_ioc_constraints(feature_win, feature_full, indToUse_window(ind_windowCount, 1), param, 'doc_rmse');

        for ind_pivot = 1:length(ccost_win)
            % regenerating the traj with the recovered weights
            currRMSE_q = 1e6;
            currRMSE_dq = 1e6;
%             q_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
%             dq_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
            
            fprintf('Solving the direct problem using pivot weights recovered (act: %u/%u)..\n', ind_pivot, length(ccost_win));
            
            if isempty(find(param.cost_functions_ioc == ind_pivot, 1))
                currRMSE_set(ind_pivot).array = 1e6;
                feature_recon_local{ind_windowCount}{ind_pivot} = [];
            else 
                try
                    if ccost_win(ind_pivot) ~= 0
                        c_pivot = ccost_inverse_pivot{ind_pivot}*ccost_win(ind_pivot); % account for non-zero pivots in the real case
                    else
                        c_pivot = ccost_inverse_pivot{ind_pivot};
                    end
                    
                    if runFirstDOCOnly && ind_pivot > 1
                        paramLa = update_intermed_ind(param, param.doc_pivot.knot_count);
                        y_opt_mat = zeros(size(paramLa.intermed_ind)) + 1e6;
                        splineFit.q = calc_direct_spline(y_opt_mat, paramLa);
                        feature_recon_local{ind_windowCount}{ind_pivot} = calc_features(splineFit, feature_win, paramLa);
                        direct_check_flags{ind_pivot}.exitflag = -3;
                        direct_check_flags{ind_pivot}.J = 1e6;
                        direct_check_flags{ind_pivot}.c = 1e6;
                        direct_check_flags{ind_pivot}.ceq_x = 1e6;
                        direct_check_flags{ind_pivot}.ceq_dx = 1e6;
                        direct_check_flags{ind_pivot}.ceq_ddx = 1e6;
                    else
                        [feature_recon_local{ind_windowCount}{ind_pivot}, ~, ~,  direct_check_flags{ind_pivot}] = main_direct(c_pivot, cconst, feature_win, param, param.doc_pivot);
                    end
                    
                    % calculate the one with the minimal RMSE so we can 
                    % decide which pivot is correct
                    currRMSE_set(ind_pivot) = calc_rmse(rmse_fct, feature_win, feature_recon_local{ind_windowCount}{ind_pivot}, param);
                    currRMSE_q = currRMSE_set(ind_pivot).q;
                    currRMSE_dq = currRMSE_set(ind_pivot).dq;
                    
%                     q_recon{ind_windowCount}{ind_pivot} = feature_recon_local{ind_windowCount}{ind_pivot}.q;
%                     dq_recon{ind_windowCount}{ind_pivot} = feature_recon_local{ind_windowCount}{ind_pivot}.dq;
                catch err
                    errMessageTmp = regexp(err.message,',','split'); % error messages with commas in it
                    errMessageFirstComma = errMessageTmp{1};
                    errMessage = [errMessageFirstComma ' - ' err.stack(1).file];
                    fprintf('Error START: [%s]\n', errMessage);
                    
                    for err_ind = 1:length(err.stack)
                        %             print the error messages to screen
                        err.stack(err_ind)
                    end
                    
                    fprintf('Error END: [%s]\n', errMessage);
                end
            end
            
            rmse(ind_pivot) = currRMSE_q;
            rmse_dq(ind_pivot) = currRMSE_dq;
            rmse_array(:, ind_pivot) = currRMSE_set(ind_pivot).array;
            output_inverse{ind_windowCount}(ind_pivot).rmse = currRMSE_q; % saving this for postpoc later
        end
        
%         docj_note = J_optsim_array_array{ind_windowCount};
%         iocj_note = output_inverse{ind_windowCount}(1).J_out_array';
%         resNormArray_note = [output_inverse{ind_windowCount}(:).resnorm]';
%         rmse_note = rmse(:);
%         rmse_set = [currRMSE_set(:).array];
        
        [~, minRmseInd] = min(rmse);
        minRmseVal(ind_windowCount) = rmse(minRmseInd);
        
        case 'resnorm_rmse1'
            fprintf('Solving the direct problem using min(resnorm) recovered... \n');
    
            % we don't want to run the RMSE checks.
            resNormArray = [output_inverse{ind_windowCount}(:).resnorm];
%             resNormArray(resNormArray == 0) = Inf;
            [~, minRmseInd] = min(resNormArray);
            
            for ind_pivot = 1:length(ccost_inverse_pivot)
                switch ind_pivot
                    case minRmseInd
                        c_pivot = ccost_inverse_pivot{ind_pivot};
                        
                        try
                            [feature_recon_local{ind_windowCount}{ind_pivot}, ~, ~,  direct_check_flags{ind_pivot}] = main_direct(c_pivot, cconst, feature_win, param, param.doc_pivot);
                            
%                             [ccost_inverse_pivot_doublecheck, J_inverse_doublecheck, output_inverse_doublecheck{ind_windowCount}, ~, ...
%                                 J_costfull_inverse_doublecheck{ind_windowCount}, J_costconst_check_doublecheck{ind_windowCount}] = ...
%                                 main_inverse2(feature_recon_local{ind_windowCount}{ind_pivot}, ccost_win, cconst, param, param.ioc, ccost_win_init, [], saveprefix);
%                             
%                             lala = 1; % break here
                            
                        catch err
                            errMessageTmp = regexp(err.message,',','split'); % error messages with commas in it
                            errMessageFirstComma = errMessageTmp{1};
                            errMessage = [errMessageFirstComma ' - ' err.stack(1).file];
                            fprintf('Error START: [%s]\n', errMessage);
                            
                            for err_ind = 1:length(err.stack)
                                %             print the error messages to screen
                                err.stack(err_ind)
                            end
                            
                            fprintf('Error END: [%s]\n', errMessage);
                            
                            rmse(ind_pivot) = 1e6;
%                             q_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
%                             dq_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
                            feature_recon_local{ind_windowCount}{ind_pivot} = feature_win;
                            currRMSE_set(ind_pivot) = calc_rmse(rmse_fct, feature_win, feature_recon_local{ind_windowCount}{ind_pivot}, param);
                            rmse_array(:, ind_pivot) = currRMSE_set(ind_pivot).array;
                            
                            direct_check_flags{ind_pivot}.exitflag = 0;
                            direct_check_flags{ind_pivot}.J = 0;
                            direct_check_flags{ind_pivot}.c = 0;
                            direct_check_flags{ind_pivot}.ceq_x = 0;
                            direct_check_flags{ind_pivot}.ceq_dx = 0;
                            direct_check_flags{ind_pivot}.ceq_ddx = 0;
                            direct_check_flags{ind_pivot}.J_breakdown.report = 0;
                        end
                            
                            
                        currRMSE_set(ind_pivot) = calc_rmse(rmse_fct, feature_win, feature_recon_local{ind_windowCount}{ind_pivot}, param);
                        currRMSE_q = currRMSE_set(ind_pivot).q;
                        currRMSE_dq = currRMSE_set(ind_pivot).dq;
                        
%                         q_recon{ind_windowCount}{ind_pivot} = feature_recon_local{ind_windowCount}{ind_pivot}.q;
%                         dq_recon{ind_windowCount}{ind_pivot} = feature_recon_local{ind_windowCount}{ind_pivot}.dq;
%                         ddq_recon{ind_windowCount}{ind_pivot} = feature_recon_local{ind_windowCount}{ind_pivot}.ddq;
                         
                        rmse(ind_pivot) = currRMSE_q;
                        minRmseVal(ind_windowCount) = currRMSE_q;
                        rmse_array(:, ind_pivot) = currRMSE_set(ind_pivot).array;
                           
                    otherwise
                        rmse(ind_pivot) = 1e6;
%                         q_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
%                         dq_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
                        feature_recon_local{ind_windowCount}{ind_pivot} = feature_win;
                        currRMSE_set(ind_pivot) = calc_rmse(rmse_fct, feature_win, feature_recon_local{ind_windowCount}{ind_pivot}, param);
                        rmse_array(:, ind_pivot) = currRMSE_set(ind_pivot).array;
                        
                        direct_check_flags{ind_pivot}.exitflag = 0;
                        direct_check_flags{ind_pivot}.J = 0;
                        direct_check_flags{ind_pivot}.c = 0;
                        direct_check_flags{ind_pivot}.ceq_x = 0;
                        direct_check_flags{ind_pivot}.ceq_dx = 0;
                        direct_check_flags{ind_pivot}.ceq_ddx = 0;
                        direct_check_flags{ind_pivot}.J_breakdown.report = 0;
                end
                
                output_inverse{ind_windowCount}(ind_pivot).rmse = rmse(ind_pivot);
            end
            
        case 'resnorm'
            % use resnorm to determine the pivot, then determine how good
            % of a fit that pivot it using DOC
                % we don't want to run the RMSE checks.
%             resNormArray = [output_inverse{ind_windowCount}(:).resnorm];
            resNormArray = [output_inverse{ind_windowCount}(:).resnorm_lsqlin];
%             resNormArray(resNormArray == 0) = Inf;
            [~, minRmseInd] = min(resNormArray);
            
            for ind_pivot = 1:length(ccost_inverse_pivot)
                rmse(ind_pivot) = 0;
%                 q_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
%                 dq_recon{ind_windowCount}{ind_pivot} = zeros(size(feature_win.q));
                feature_recon_local{ind_windowCount}{ind_pivot} = feature_win;
%                 currRMSE_set(ind_pivot) = calc_rmse(rmse_fct, feature_win, feature_recon_local{ind_windowCount}{ind_pivot}, param);
%                 rmse_array(:, ind_pivot) = currRMSE_set(ind_pivot).array;
                
                currRMSE_set(minRmseInd).q = 0;
                rmse_array(:, ind_pivot) = 0;
                direct_check_flags{ind_pivot}.exitflag = 0;
                direct_check_flags{ind_pivot}.J = 0;
                direct_check_flags{ind_pivot}.c = 0;
                direct_check_flags{ind_pivot}.ceq_x = 0;
                direct_check_flags{ind_pivot}.ceq_dx = 0;
                direct_check_flags{ind_pivot}.ceq_ddx = 0;
                direct_check_flags{ind_pivot}.J_breakdown.report = 0;
            end
            
            minRmseVal(ind_windowCount) = 0;
    end
    
    [minVal, minInd] = min(rmse_array, [], 2);
    minIndCount= [];
    for ind_minInd = 1:length(ccost_win)
        minIndCount(ind_minInd) = sum(minInd == ind_minInd);
    end
    
    lambda_array = [output_inverse{ind_windowCount}(:).lambda];

    % save the one with the smallest RMSE, and other corresponding data
    t_recon_plot_array{ind_windowCount} = feature_win.t;
%     indknot_recon_plot_array{ind_windowCount} = param.intermed_ind_set;
%     indqconst_recon_plot_array{ind_windowCount} = param.const_x;
%     inddqconst_recon_plot_array{ind_windowCount} = param.const_dx;
    q_recon_plot_array{ind_windowCount} = feature_recon_local{ind_windowCount}{minRmseInd}.q;
    dq_recon_plot_array{ind_windowCount} = feature_recon_local{ind_windowCount}{minRmseInd}.dq;
    ddq_recon_plot_array{ind_windowCount} = feature_recon_local{ind_windowCount}{minRmseInd}.ddq;
    
    minRmseIndArray(ind_windowCount) = minRmseInd;
    
    switch param.segment_only_windows
        case 'none'
            frontEdgeInd = indToUse_window(ind_windowCount, 2)-param.win_shift+1;
            if frontEdgeInd < 1
                frontEdgeInd = 1;
            end
            
        otherwise
            frontEdgeInd = indToUse_window(ind_windowCount, 1);
    end
    
    ind_toUseT = frontEdgeInd:indToUse_window(ind_windowCount, 2); % the time length covered is the same as the window inserted
    lenToUseT = length(ind_toUseT);
    t_rmse = [t_rmse feature_full.t(:, ind_toUseT)]; % cover the increment between the previous window and this one
    dofsUsedInIOC = zeros(size(ccost_win));
    dofsUsedInIOC(param.cost_functions_ioc) = param.cost_functions_ioc;
    
    ratioAll = horzcat(J_inverse{:});
    ratioRMSE{ind_windowCount} = J_inverse{minRmseInd};
    
    ratioAll_plot =   [ratioAll_plot;   repmat(ratioAll,                           [lenToUseT 1])];
    ratioRMSE_plot =  [ratioRMSE_plot;  repmat(ratioRMSE{ind_windowCount},         [lenToUseT 1])];
    rmseMin_plot =    [rmseMin_plot;    repmat(rmse(minRmseInd),                   [lenToUseT 1])];
    condMin_plot =    [condMin_plot;    repmat(output_inverse{ind_windowCount}(minRmseInd).condA,   [lenToUseT 1])];
    resnormMin_plot = [resnormMin_plot; repmat(output_inverse{ind_windowCount}(minRmseInd).resnorm, [lenToUseT 1])];
%     q_jh_range_plot = [q_jh_range_plot; repmat(output_inverse{ind_windowCount}(minRmseInd).corrJHqrange_overall, [lenToUseT 1])];
%     q_jh_absmax_plot =[q_jh_absmax_plot;repmat(output_inverse{ind_windowCount}(minRmseInd).corrJHabsmax_overall, [lenToUseT 1])];
    
    dofsUsedInIOC_plot = [dofsUsedInIOC_plot; repmat(dofsUsedInIOC,     [lenToUseT 1])];
    minRmseInd_plot =    [minRmseInd_plot; repmat(minRmseInd,           [lenToUseT 1])]; 
    
    rmseAll_plot =    [rmseAll_plot;    repmat(rmse,                                         [lenToUseT 1])];
    condAll_plot =    [condAll_plot;    repmat([output_inverse{ind_windowCount}(:).condA],   [lenToUseT 1])];
    resnormAll_plot = [resnormAll_plot; repmat([output_inverse{ind_windowCount}(:).resnorm], [lenToUseT 1])];
    resnormAll_lsqlin_const_plot =   [resnormAll_lsqlin_const_plot;   repmat([output_inverse{ind_windowCount}(:).resnorm_lsqlin], [lenToUseT 1])];
%     resnormAll_lsqlin_unconst_plot = [resnormAll_lsqlin_unconst_plot; repmat([output_inverse{ind_windowCount}(:).resnorm_lsqlin_unconst], [lenToUseT 1])];
    resnormAll_lsqlin_const_RMSE_plot =   [resnormAll_lsqlin_const_RMSE_plot;   repmat([output_inverse{ind_windowCount}(minRmseInd).resnorm_lsqlin], [lenToUseT 1])];
%     resnormAll_lsqlin_unconst_RMSE_plot = [resnormAll_lsqlin_unconst_RMSE_plot; repmat([output_inverse{ind_windowCount}(minRmseInd).resnorm_lsqlin_unconst], [lenToUseT 1])];
    
    c_ioc_constrain_plot = [c_ioc_constrain_plot; repmat(output_inverse{ind_windowCount}(minRmseInd).c_recovered, [lenToUseT 1])];
%     c_ioc_unconstrain_plot = [c_ioc_unconstrain_plot; repmat(output_inverse{ind_windowCount}(minRmseInd).c_out_pivot_unconst, [lenToUseT 1])];
    
    for ii_pivotInspect = 1:length(ccost_inverse_pivot)
        ccost_offset = ccost_win(ii_pivotInspect);
        if ccost_offset == 0
            ccost_offset = 1;
        end
        
        weight_offset_val{ii_pivotInspect} = ccost_inverse_pivot{ii_pivotInspect}*ccost_offset;
    end
    
    ccostRMSE_select = weight_offset_val{minRmseInd};
    ccostRMSE_plot = [ccostRMSE_plot;   repmat(ccostRMSE_select,               [lenToUseT 1])];
    ccostAll_plot = [ccostAll_plot;     repmat(horzcat(weight_offset_val{:}),  [lenToUseT 1])];

    if exist('output_optsim', 'var')
        doc_optsim_win = output_optsim(output_optsim_mapping_full(indToUse_window(ind_windowCount, 1)));
    else
        doc_optsim_win.grad = 0;
        doc_optsim_win.exitflag = 0;
        doc_optsim_win.J_array = 0;
        doc_optsim_win.J_contrib = zeros(size(cost_function_names));
        doc_optsim_win.c = 0;
        doc_optsim_win.J = 0;
%         doc_optsim_win.J_breakdown = 0;
        doc_optsim_win.ceq_x = 0;
        doc_optsim_win.ceq_dx = 0;
        doc_optsim_win.ceq_ddx = 0;
%         doc_optsim_win.J_breakdown.report = zeros(size(direct_check_flags{1}.J_breakdown.report));
    end
    
    report2; format long % C_report
    ds = cell2dataset(C_report); 
    targetInstanceCsv = fullfile(outputPath_intermedCsv, ...
        [currInstName '_log_' num2str(ind_windowCount) '_' num2str(feature_win.t(1), '%1.1f') '_'  num2str(feature_win.t(end), '%1.1f') '.csv']);
    export(ds,'file',targetInstanceCsv,'delimiter',',');
    
    % write the correlation matrix to the file
%     output_inverse(ind_pivotset).corrValFull
    fileID = fopen(targetInstanceCsv,'a');
    
    writeMatrixToFile(fileID, '\n J_costconst = \n', J_costconst_check{ind_windowCount});
    writeMatrixToFile(fileID, '\n mean(abs(J cost inv full)) = \n', mean(abs(J_costfull_inverse{ind_windowCount})));
%     writeMatrixToFile(fileID, 'corrcoef(J_costconst)', corrcoef(J_costconst_check{ind_windowCount}));
% writeMatrixToFile(fileID, '\n correlation matrix (all) = \n', output_inverse{ind_windowCount}(1).corrValMatrix_All);
% writeMatrixToFile(fileID, '\n\n correlation matrix (nan) = \n', output_inverse{ind_windowCount}(1).corrValMatrix_NaN);
% writeMatrixToFile(fileID, '\n J cost inv full correlation matrix = \n', J_costfull_inverse{ind_windowCount});

    fclose(fileID);
    
    [motionEnd, motionWhole, ~, ~] = ...
        checkStandingVSMovement(feature_full, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2), ...
        segmentInfo, ratioRMSE{ind_windowCount}, [], []);
    
    writeToInstanceCSV(fullfile(outputPath, ['' currInstName '_tally.csv']), cost_function_names, indToUse_window(ind_windowCount, 1), indToUse_window(ind_windowCount, 2), ...
        set_ioc_check_struct, ccost_win, doc_optsim_win, ccost_inverse_pivot, ...
        J_inverse, currRMSE_set, output_inverse, ind_windowCount, minRmseInd, direct_check_flags, motionWhole, param);
    
    if runSettings.saveResults && mod(ind_windowCount, safetySave) == 0
        intermedMatPath{end+1} = fullfile(outputPath_intermedMat, [currInstName '_' num2str(ind_windowCount)]);
        save(intermedMatPath{end});
    end
end

% if runSettings.saveResults 
%     intermedMatPath{end+1} = fullfile(outputPath_intermedMat, [currInstName '_' num2str(ind_windowCount)]);
%     save(intermedMatPath{end});
% end

rmse_report.t = feature_full.t;
rmse_report.maxminRMSE = max(minRmseVal);
rmse_report.meanRMSE = mean(minRmseVal);
rmse_report.rangeRMSE = range(minRmseVal);
rmse_report.stdRMSE = std(minRmseVal);

elapsedTime = toc;
toc

try
    if runSettings.saveResults
        fprintf('Saving output at %s\n', fullfile(outputPath, [currInstName '_end']));
        param = rmfield(param, 'model');
        
		saveVarGenerate();
        save(fullfile(outputPath, [currInstName '_end']), 'saveVar');
        
        % once this is done, remove all the intermediate files
        for ind_intermed = 1:length(intermedMatPath)
            delete([intermedMatPath{ind_intermed} '.mat']);
        end
    end
catch err
    errMessageTmp = regexp(err.message,',','split'); % error messages with commas in it
    errMessageFirstComma = errMessageTmp{1};
    errMessage = [errMessageFirstComma ' - ' err.stack(1).file];
    fprintf('Error START: [%s]\n', errMessage);
    
    for err_ind = 1:length(err.stack)
        %             print the error messages to screen
        err.stack(err_ind)
    end
    
    fprintf('Error END: [%s]\n', errMessage);
end


beep
beep
beep
beep
beep

end