function output = zeroScaleSignal(input, zeroScale)
    % scale it so it only outputs values between 0 and 1
    
    output = zeros(size(input));
    if zeroScale
        % offset to [0, n]
        for i = 1:size(input, 2)
            output(:, i) = input(:, i)  - min(input(:, i));
            
            output(:, i) = output(:, i) / max(output(:, i));
        end
        
        % keep the relative scaling between the DOFs
%         output = output / max(max(output));
    else
        output = input;
    end