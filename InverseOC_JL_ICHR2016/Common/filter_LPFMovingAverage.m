function lpfVal = filter_LPFMovingAverage( prefilterVal, filterBracket )
    % PURPOSE - Causal LPF, implemented via evenly-weighted averaging of
    % 'filterBracket' number of previous data. If previous data does not 
    % exist (ie the first value passed in), averaging is done with the 
    % most number of items possible
    
    
    array_length = length(prefilterVal);
    lpfVal = zeros(size(prefilterVal));
    
    for counter = 1:array_length
        clicker = filterBracket;
        if (counter - filterBracket) < 1
            % So if we take 'filterBracket' number of previous data, we
            % would attempt to read from negative indices of the array...
            clicker = counter - 1;
        end
        
        lpfVal(counter, :) = mean(prefilterVal(counter - clicker:counter, :), 1);
    end
end
