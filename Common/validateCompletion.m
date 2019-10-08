function [rank, completed, errorCode] = validateCompletion(Hhat, gamma, delta, dimWeights)
% Evaluates recovery matrix rank to determine whether it is possible to
% recover weights or more points are necessary

    s = svd(Hhat);
    s = sort(s);
    
    % These are the conditions outlined in paper. However, Wanxin's
    % implementation uses slightly different conditions
    rank = abs(s(2)/s(1));
    
    errorCode(1) = s(1) < delta;
    errorCode(2) = rank > gamma;  
    errorCode(3) = hPassFct1(Hhat, dimWeights);
%         if (s(1)<delta &&  rank > gamma && rows >= dimWeights * cols)
    if (sum(errorCode) == length(errorCode)) % passed all the tests
        completed = 1;
    else
        completed = 0;
    end
end

