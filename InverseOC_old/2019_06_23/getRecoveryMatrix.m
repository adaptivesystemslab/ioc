function [H1, H2]= getRecoveryMatrix(iocInstance, prevH1, prevH2, x0, u0, x, u, dt)
    
    numStates = size(x,2);
    dofs = numStates/2;
   
    % Get features derivatives
    [~, ~, px, pu] = iocInstance.getDerivativesNewObservation(u0, x(:,1:dofs),...
        x(:,dofs+1:end), dt); 
    dp_du = pu';
    dp_dx = px';
    % Get dynamic derivatives wrt to control
    [~, fu, ~, ~] = iocInstance.getDerivativesNewObservation(u0, x0(:,1:dofs),...
        x0(:,dofs+1:end), dt); % I'm supposing that in x vector, 1-4 cols are angles and 5-8 are velocities
    df_du = fu';
    % Get dynamic derivatives wrt to state
    % Derivatives x1
    [fx, ~, ~, ~] = iocInstance.getDerivativesNewObservation(u, x(:,1:dofs),...
        x(:,dofs+1:end), dt);
    df_dx = fx';
    
    
    if (isempty(prevH1) && isempty(prevH2))
        H1=df_du*dp_dx+dp_du;
        H2=df_du*df_dx;
    else
        H1 = [prevH1+prevH2*dp_dx;
              df_du*dp_dx+dp_du];
        H2 = [prevH2*df_dx;
              df_du*df_dx];
    end
end
