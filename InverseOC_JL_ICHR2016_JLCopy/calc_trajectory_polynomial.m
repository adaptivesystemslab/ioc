function [q, dq, polyParam] = calc_trajectory_polynomial(t, y)
    Nh = floor(length(t)/2);
    Nf = length(t);
    
    for j = 1:size(y, 1)
        y_use = y(j, :);
        
        p = polyfit(t, y_use, 5);
        q(j, :) = polyval(p, t);
        dq(j, :) = calcDeriv(q(j, :), mean(diff(t)));
        
        polyParam(j, :) = p;
    end
end