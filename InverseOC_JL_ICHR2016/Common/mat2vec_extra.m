function [state_vector, feature_count, extra_array_length] = mat2vec_extra(state_array, extra_array)
    % this function converts a matrix form 'state_array' into a vector form
    % 'state_vector'. it is designed to be used in conjunction with
    % 'vec2mat' to allow for easy conversion between the two forms
    % if an 'extra_array' is passed in, it is concatentated to the
    % state_vector
    
    extra_array = extra_array(:);
    extra_array_length = length(extra_array);
    
    [state_vector, feature_count] = mat2vec(state_array);
    state_vector = [extra_array; state_vector];
end