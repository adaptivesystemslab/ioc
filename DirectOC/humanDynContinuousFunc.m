function phaseout = humanDynContinuousFunc(input, weights, iocObject)
    
    w = weights'; w = w/sum(w);
    
    x = input.phase.state; 
    u = input.phase.control; 
    t = input.phase.time;
    
    dofs = size(u,2);
    
    q = x(:,1:dofs);
    dq = x(:,dofs+1:end);
           
    for i = 1:size(q,1)     
        ddq(i,:) = iocObject.dynamicModel.forwardDynamicsQDqTau(q(i,:), dq(i,:),u(i,:));
    end
    
    phaseout.dynamics = [dq, ddq];
 
    
    %Compute features
%     normFeatures = iocObject.calcFeatures(x, u)/1;
%      normFeatures = features ./ max(features(:));
%     phaseout.integrand = sum(w.*normFeatures,2);
     phaseout.integrand=100*u(:,1).^2+1*u(:,2).^2;
end