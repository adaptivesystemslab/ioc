function error = computeError(estimatedCost, trueCost)
%COMPUTEERROR Summary of this function goes here
%   Detailed explanation goes here
    l=size(estimatedCost,2);
    error=[];
    estimatedCost = estimatedCost'; % Added this line since I only have a single vector
    for i=1:l
        est=estimatedCost(i,:);
        c=dot(est,trueCost)/dot(est,est);
        e=norm(c*est-trueCost)/norm(trueCost);
        error(end+1)=e;
    end
end

