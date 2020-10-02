function output = bw_dualfilt(filterVal, filterFreq, filterOrder, type)
    % implements a double pass butterworth filter
    [b, a] = butter(filterOrder,filterFreq,type);

    for x = 1:size(filterVal, 2)
        output(:, x) = filtfilt(b, a, filterVal(:, x));
    end
end