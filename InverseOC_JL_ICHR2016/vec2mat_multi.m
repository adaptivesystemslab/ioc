function [state_array, extra_vector] = vec2mat(state_vector, m2v_struct)
    % this function converts a vector form 'state_vector' into a matrix,
    % given that we know what the statecount (the width of the matrix) a
    % priori. it is designed to be used in conjunction with vec2mat to 
    % allow for easy conversion between the two forms
    
    feature_count = m2v_struct.feature_count;
    feature_width = m2v_struct.feature_width;
    entry_count = m2v_struct.entry_count;
    
    if ~isfield(m2v_struct, 'extra_array_length')
        extra_vector = [];
    else
        extra_array_length = m2v_struct.extra_array_length;
        
        extra_vector = state_vector(1:extra_array_length);
        state_vector = state_vector(extra_array_length+1:end);
    end
    
    state_array_temp = reshape(state_vector, feature_count, size(state_vector, 1)/feature_count); % size(2), size(1)
    for ind = 1:entry_count
        currIndX = (1:feature_count);
        currIndY = (ind-1)*feature_width + (1:feature_width);
        state_array{ind} = state_array_temp(currIndX, currIndY);
    end
end

