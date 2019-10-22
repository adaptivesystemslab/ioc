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
    features = iocObject.calcFeatures(x, u);
    normFeatures = features ./ max(features(:));
    phaseout.integrand=normFeatures*w';
    
%     figure(1)
%     clf
%     plot(t,q);
%     drawnow
%     pause;
end