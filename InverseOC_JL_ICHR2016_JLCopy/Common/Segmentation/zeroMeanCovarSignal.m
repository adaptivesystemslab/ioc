function [output, inputMean, inputStd] = zeroMeanCovarSignal(input, zeroMean, zeroCovar, meanToOffset, stdDevToOffset)
    % shift and scale the signal so that it has 0 mean and 1 covar, or
    % scale it so it only outputs values between 0 and 1
    
    if ~exist('zeroMean', 'var')
        zeroMean = 1;
    end
    
    if ~exist('zeroCovar', 'var')
        zeroCovar = 0;
    end
    
    if zeroMean
        if exist('meanToOffset', 'var') || ~isempty(meanToOffset)
            inputMean = meanToOffset;
        else
            inputMean = mean(input);
        end
    else
        inputMean = zeros(1, size(input, 2));
    end
    
    if zeroCovar
        if exist('stdDevToOffset', 'var') || ~isempty(stdDevToOffset)
            inputStd = stdDevToOffset;
        else
            inputStd = std(input);
        end
    else
        inputStd = ones(size(inputMean));
    end
    
    if zeroMean || zeroCovar
        output = zeros(size(input));
        for i = 1:size(input, 2)
            output(:, i) = (input(:, i) - inputMean(i)) / inputStd(i);
        end
    else
        % is not actually doing any scaling since none of the flags are set
        output = input;
    end