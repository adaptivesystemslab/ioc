function [J, J_contrib, J_array, J_debug] = calc_direct_cost(c_cost, feature_use, param)

    if isempty(c_cost)
        c_cost = ones(size(param.cost_function_names))';
    end
    
    % feature normalization only needs to occur before cost calculation
%     feature_normalized = normalize_features(feature_use, param);
    feature_normalized = feature_use;
    
    % calculate the individual cost function values, based on the name 
    % array in main.m
    J_array = zeros(size(param.cost_function_names));
    for i = 1:length(param.cost_function_names)
        J_array(i) = calc_cost_function(param.cost_function_names{i}, feature_normalized, param);
    end
    
    if 0
        ind_set{1} = 1:length(param.t_spline);
        ind_set{2} = 1:floor(length(param.t_spline)/2);
        ind_set{3} = 1:floor(length(param.t_spline)/2)+1;
        ind_set{4} =  floor(length(param.t_spline)/2)+1:length(param.t_spline);
        
        switch length(param.intermed_ind)
            case 5
                ind_set{5} = param.intermed_ind(2):param.intermed_ind(4);
        
            otherwise
                ind_set{5} = 1:1;
        end
        
        for j = 1:length(ind_set)
            for i = 1:length(param.cost_function_names)
                J_array_set{j}(i, 1) = calc_cost_function_ind(param.cost_function_names{i}, feature_normalized, param, ind_set{j});
            end
            
            J_debug_ind(:, j) = [ind_set{j}(1); ind_set{j}(end);];
        end
        
        J_debug.ind = J_debug_ind;
        J_debug.J_agg = [J_array_set{:}];
        
        la = [J_debug.ind' J_debug.J_agg']';
        J_debug.report = la(:);
    else
        J_debug.ind_all = 0;
        J_debug.ind_start = 0;
        J_debug.ind_end = 0;
        J_debug.J_array_all = 0;
        J_debug.J_array_start = 0;
        J_debug.J_array_end = 0;
        J_debug.ind = [0 0 0 0 0 0];
        J_debug.J_agg = [];
        J_debug.la = [];
        J_debug.report = [];
    end
    
    % and apply the cost weight
    J_cost = c_cost .* J_array;
    J = sum(J_cost);
    
    % calculate the contribution ratio
    J_contrib = J_cost / J * 100;
end