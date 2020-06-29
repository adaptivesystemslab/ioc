function [weights, residual] = computeWeights(Hhat, numWeights)
% Obtain weights from recovery matrix using cvx
    sumConstraint = 1;
    [~, c] = size(Hhat);

    %cvx specification 
    cvx_begin quiet 
        variable x(c)    
        minimize (norm(Hhat*x,2))
        subject to    
        sum(x(1:numWeights))==sumConstraint; %#ok<*EQEFF> 
    cvx_end 

    weights=x(1:numWeights);
    if size(weights,1)==1
        weights=weights'; 
    end
    residual=cvx_optval/norm(x);
end
