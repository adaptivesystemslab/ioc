function [output, zf] = filter_butterworth_function(filterVal, type, pass, filterFreq, order, zi, padLength)
    % PURPOSE - provides a single-pass Butterworth low-pass filter
    % Accepts 1-dim or 2-dim (vertical vectors) data
    
    % filterVal     array to filter
    % type          filter type: low, high, band
    % pass          how many passes: single, double
    % cutoff        cutoff frequency
    % sampleRate    sample rate of the data
    % order         order of Butterworth.
    % zi            initial values for the filter
    % padLength     if non-zero, denotes the amount of pre-padding that
    %               should be added to the data to reduce rise time issues
    
    % output        filtered data
    % zf            final condition of the butterworth filter
    
    % ex: output = dualpassBW(filterVal, 30, 300, 10)

    dataDof = size(filterVal, 2);
    
    if ~exist('type', 'var')  || isempty(type)
        type = 'low';
        pass = 'doublePass';
    end
    
    if ~exist('filterFreq', 'var') || isempty(filterFreq)
        % if the function is called without filtering freq...
        filterFreq = 0.04;
        order = 2;
    end
    
    if ~exist('zi', 'var') || isempty(zi)
        zi = zeros(order, dataDof);
    end
    
    if ~exist('padLength', 'var')
        % pre-padding length, if need to match initial value
        padLength = 0;
    end

    % generate the filter parameters
    [b, a] = butter(order,filterFreq,type);

    % generate temp variable to pass out
    output = zeros(size(filterVal));
    zf = zeros(order, dataDof);
    
    for x = 1:dataDof
        arrayToFilter = filterVal(:, x);
        
        if padLength
            % insert a padding to 'wind up' the filter           
            arrayToFilter = padarray(arrayToFilter, padLength, arrayToFilter(1), 'pre');
        end
        
        [dofOutput, dofZf] = filterFunction(arrayToFilter, pass, b, a, zi(:, x));
        
        if padLength
            % then crop out the 'wind up' rise time
            dofOutput = dofOutput(padLength+1:end);
        end
        
        output(:, x) = dofOutput;
        
        if sum(zi) > 0
            zf(:, x) = dofZf;
        else
            zf = [];
        end
    end
    
%     % testing
%     plot(filterVal(:, 1), 'b');
%     hold on
%     plot(dofOutput(:, 1), 'r');
       
end

function [output, zf]  = filterFunction(filterVal, pass, b, a, zi)

    if sum(zi) > 0
        if strcmpi(pass, 'singlePass')
            [output, zf] = filter(b,a,filterVal, zi);
        elseif strcmpi(pass, 'doublePass')
            [output, zf]  = filtfilt(b,a,filterVal, zi);
        end   
    else
        if strcmpi(pass, 'singlePass')
            [output] = filter(b,a,filterVal);
        elseif strcmpi(pass, 'doublePass')
            [output]  = filtfilt(b,a,filterVal);
        end   
        zf = [];
    end
end