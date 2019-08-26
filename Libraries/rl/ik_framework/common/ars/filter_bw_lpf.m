function output = filter_bw_lpf(filterVal, cutoff, order)
    % PURPOSE - provides a single-pass Butterworth low-pass filter
    % Accepts 1-dim or 2-dim (vertical vectors) data
    
    % filterVal     array to filter
    % cutoff        cutoff frequency
    % order         order of Butterworth. Dual-pass system, so the order 
    %               is effectively double of this
    
    % ex: output = dualpassBW(filterVal, 30, 300, 10)

    if ~exist('cutoff', 'var')
        % if the function is called without filtering freq...
        cutoff = 0.04;
        order = 5;
    end
    
    [b, a] = butter(order,cutoff,'low');
    output = filtfilt(b,a,filterVal);
    
%     if max(size(filterVal)) == 1
%         % if the filterVal matrix is only one dimensional
%         forwardFilter = filterFunction(filterVal, filterFreq, 'low', order);
%         output = forwardFilter;
%     else    
%         output = zeros(size(filterVal));
%         for x = 1:size(filterVal, 2)
%             forwardFilter = filterFunction(filterVal(:, x), filterFreq, 'low', order);
%             output(:, x) = forwardFilter;
%         end
%     end
end