% load a specific .mat file and using the fixed window cf functions,
% recover the ind window values
matPath = 'C:\Users\jf2lin\Downloads\TransferXL-00j4BNd2YW3mTp\Subject01_SingleArm60Time_4071_7071_end.mat';
if ~exist('saveVar', 'var')
    load(matPath);
end

addpath('Common');

lenWin = length(saveVar.feature_win_save);

% and feat, because Jc = sum(sum(feat .^ 2)) / len
for ind_win = 1:lenWin
   loadStuff(saveVar, ind_win);
end

% J = sum w * J/J_norm, and give you each of the individual values


function loadStuff(saveVar, ind_win)
    % the windowed data from the original trajectory (NOT POST SPLINE)
    feature_opt = saveVar.feature_win_save{ind_win}; 
    c_use = ones(size(saveVar.cost_function_names, 2), 1);
    
    param = saveVar.param;
    param.cost_function_names = saveVar.cost_function_names;
%     param.intermed_ind = [1 21 41 61]; 
%     param.spline_order = 5;
%     param.splineType = saveVar.runSettings.variableFactors.splineType;
%     param.splineSlave = saveVar.runSettings.variableFactors.splineSlave;
%     param.dofsFull = size(feature
    
    
%     % pull out knots from the traj, then generate a spline-based traj from that
%     y_opt_base{1} = feature_opt.q(:, param.intermed_ind);
%     y_opt_base{2} = feature_opt.dq(:, param.intermed_ind);
%     y_opt_base{3} = feature_opt.ddq(:, param.intermed_ind);
%     [y_opt_vec, param.matVec_struct] = mat2vec_multi(y_opt_base); 
%     y_opt_mat = vec2mat_multi(y_opt_vec, param.matVec_struct);
%     curr_spline_fit.q =  calc_direct_spline(y_opt_mat, param);
%     
%     % using the spline, calculate the features
%     feature_use = calc_features(curr_spline_fit, [], param);
    
    % compute the cost function values corresponding to this combination.
    % that is, how does the cost function change when we change the
    [J, J_contrib, J_array, J_cost_array] = calc_direct_cost(c_use, feature_opt, param);
end