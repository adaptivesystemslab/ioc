function [J, J_contrib, J_array] = calc_direct_cost_syms(c_cost, feature_use, param)

    if isempty(c_cost)
        c_cost = ones(size(param.cost_function_names));
    end
    
    % feature normalization only needs to occur before cost calculation
    feature_normalized = normalize_features(feature_use, param);

    % calculate the individual cost function values, based on the name 
    % array in main.m
    J_array = {};
    for i = 1:length(param.cost_function_names)
        J_array{i} = calc_cost_function_syms(param.cost_function_names{i}, feature_normalized, param);
    end
    
    J = sym(0);
    for i = 1:length(param.cost_function_names)
        J = J + c_cost(i) * J_array{i};
    end
    
    % calculate the contribution ratio
    J_contrib = 0;
end

