function [H1, H2]= getRecoveryMatrix(iocInstance, prevH1, prevH2, x, u, dt)    
    % Get dynamic derivatives wrt to control
    % Get dynamic derivatives wrt to state
    [fx, fu, px, pu] = iocInstance.getDerivativesNewObservation(x, u);
    df_dx = fx';
    df_du = fu';
    dp_dx = px';
    dp_du = pu';
    
    [H1, H2] = assembleH1H2(df_dx, df_du, dp_dx, dp_du, prevH1,prevH2);
end