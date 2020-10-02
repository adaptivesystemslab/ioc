function [state_vector, m2v_struct] = mat2vec_multi_syms(state_array_in, extra_array)
    % this function converts a matrix form 'state_array' into a vector form
    % 'state_vector'. it is designed to be used in conjunction with
    % 'vec2mat' to allow for easy conversion between the two forms
    
    % if state array is a cell, assume that all state array is the same
    % size as the first entry
    
    if nargin == 1
        extra_array = [];
    end
    
    if ~iscell(state_array_in)
        state_array{1} = state_array_in;
    else
        state_array = state_array_in;
    end
        
    feature_count = size(state_array{1}, 1);
    feature_width = size(state_array{1}, 2);
    entry_count = length(state_array);
    extra_array_length = length(extra_array);
    
    feature_space = feature_width*feature_count;
    state_vector = sym(zeros(feature_space*entry_count, 1));
    for ind = 1:entry_count
        currInd = (ind-1)*feature_space + (1:feature_space);
        state_vector(currInd) = reshape(state_array{ind}, feature_space, 1);
    end
    
    state_vector = [extra_array; state_vector];
    
    m2v_struct.feature_count = feature_count;
    m2v_struct.feature_width = feature_width;
    m2v_struct.entry_count = entry_count;
    m2v_struct.extra_array_length = extra_array_length;
end