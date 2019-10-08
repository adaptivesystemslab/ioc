function phaseout = parameterizedDynamics(input, weights, iocObject)
    w = weights'; w = w/sum(w);
    
    x = input.phase.state; 
    u = input.phase.control; 
    t = input.phase.time;
    
    % Compute current time increase and update value in iocInstance
    dt = t(2) - t(1);
    iocObject.dt = dt;
    
    dofs = size(u,2);
    
    q = x(:,1:dofs);
    dq = x(:,dofs+1:end);
           
    for i = 1:size(q,1)     
        % Set arm state
        iocObject.dynamicModel.updateState(q(i,:), dq(i,:));
        ddq(i,:) = iocObject.dynamicModel.forwardDynamics(u(i,:));
    end
    
    phaseout.dynamics = [dq, ddq];
    
    %Check bounds on states using auxiliary variables
    stateLower = input.auxdata.state.lower;
    stateUpper = input.auxdata.state.upper;
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            if ~(x(i,j) >= stateLower(j) && x(i,j) <= stateUpper(j))
                fprintf("State variable %i is out of bounds", j);
            end
        end
    end
    
    %Check bounds on control using auxiliary variables
    controlLower = input.auxdata.control.lower;
    controlUpper = input.auxdata.control.upper;
    for i = 1:size(u,1)
        for j = 1:size(u,2)
            if ~(u(i,j) >= controlLower(j) && u(i,j) <= controlUpper(j))
                fprintf("Control variable %i is out of bounds", j);
            end
        end
    end
                    
    %Compute features
    features = iocObject.calcFeatures(x, u);
    normFeatures = features ./ max(features(:));
    phaseout.integrand = sum(w.*normFeatures,2);
    
    %Check bounds on integrand
    for i = 1:size(u,1)
        if ~(phaseout.integrand(i) >= input.auxdata.integral.lower && ...
                phaseout.integrand(i) <= input.auxdata.integral.upper)
            fprintf("Integran out of bounds\n");
        end
    end
end