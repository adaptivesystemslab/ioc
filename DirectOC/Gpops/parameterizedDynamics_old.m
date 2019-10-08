function phaseout = parameterizedDynamics(input, weights, model)
    w = weights'; w = w/sum(w);
    
    x = input.phase.state; 
    u = input.phase.control; 
    
    dofs = size(u,2);
    
    q = x(:,1:dofs);
    dq = x(:,dofs+1:end);
           
    for i = 1:size(q,1)     
        % Set arm state
        model.updateState(q(i,:), dq(i,:));
        ddq(i,:) = model.forwardDynamics(u(i,:));
    end
    
    phaseout.dynamics = [dq, ddq]; 
    phaseout.integrand = sum(w.*(u.*u),2);
    %phaseout.integrand = w(1)*u1.^2 + w(2)*u2.^2 + w(3)*u3.^2 + w(4)*u4.^2;
end