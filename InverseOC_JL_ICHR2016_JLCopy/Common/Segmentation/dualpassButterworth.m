function output = dualpassButterworth(filterVal, cutoff, sampleRate, order, filtType)
    % PURPOSE - provides a dual-pass Butterworth filter via filtfilt
    % Accepts 1-dim or 2-dim (vertical vectors) data
    
    % filterVal     array to filter
    % cutoff        cutoff frequency
    % sampleRate    sample rate of the data
    % order         order of Butterworth
    % filtType      'low' 'high'
    
    % ex: output = dualpassBW(filterVal, 0.16, 0, 2, 'low')
    % ex: output = dualpassBW(filterVal, 0.90, 0, 2, 'high')

    if ~exist('cutoff', 'var')
        % if the function is called without filtering freq...
        cutoff = 0.04;
        sampleRate = 0;
        order = 5;
    end
    
    if ~exist('filtType', 'var')
        filtType = 'low';
    end
    
    if sampleRate > 0
        % if a sample rate was specified...(LEGACY SUPPORT)
        filterFreq = 2*pi*cutoff/sampleRate * (1/2);
    else
        filterFreq = cutoff;
    end
    
    if max(size(filterVal)) == 1
        % if the filterVal matrix is only one dimensional
        forwardFilter = filterFunction(filterVal, filterFreq, filtType, order);
        output = forwardFilter;
    else    
        output = zeros(size(filterVal));
        for x = 1:size(filterVal, 2)
            forwardFilter = filterFunction(filterVal(:, x), filterFreq, filtType, order);
            output(:, x) = forwardFilter;
        end
    end
end

function output = filterFunction(filterVal, filterFreq, type, filterOrder)

    [b, a] = butter(filterOrder,filterFreq,type);
    output = filtfilt(b,a,filterVal);
    
end