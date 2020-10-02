function [state_array, extra_vector_base] = vec2mat_extra(state_vector, feature_count, extra_array_length)
    % this function converts a vector form 'state_vector' into a matrix,
    % given that we know what the statecount (the width of the matrix) a
    % priori. it is designed to be used in conjunction with vec2mat to 
    % allow for easy conversion between the two forms
    
%     state_vector_base = state_vector(1:end-extra_array_length);
%     extra_vector_base = state_vector(end-extra_array_length+1:end);
    
    extra_vector_base = state_vector(1:extra_array_length);
    state_vector_base = state_vector(extra_array_length+1:end);
   
    state_array = vec2mat(state_vector_base, feature_count);
end

