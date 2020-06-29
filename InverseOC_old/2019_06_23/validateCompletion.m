function [rank, completed] = validateCompletion(Hhat, gamma, delta, dimWeights)
% Evaluates recovery matrix rank to determine whether it is possible to
% recover weights or more points are necessary

    completed = 0;
    [rows, cols] = size(Hhat);
    s = svd(Hhat);
    s = sort(s);
    
    % These are the conditions outlined in paper. However, Wanxin's
    % implementation uses slightly different conditions
    rank = abs(s(2)/s(1));
    if (s(1)<delta &&  rank > gamma && rows >= dimWeights * cols)
        completed = 1;
    end
end

