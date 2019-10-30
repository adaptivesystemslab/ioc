function dz = Dynamics(x, u, iocObject)
        
    dofs = size(u,1);
    
    x = x';
    u = u';
    
    q = x(:,1:dofs);
    dq = x(:,dofs+1:end);
           
    for i = 1:size(q,1)     
        ddq(i,:) = iocObject.dynamicModel.forwardDynamicsQDqTau(q(i,:), dq(i,:), u(i,:));
    end
    
    dz = [dq, ddq]';
    
end