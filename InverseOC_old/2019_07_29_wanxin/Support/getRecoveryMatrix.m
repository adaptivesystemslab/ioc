% function [H1, H2]= getRecoveryMatrix(iocInstance, prevH1, prevH2, x0qdq, u0, xqdq, u, dt)
%    
%     % rearrange state from [q1 q2 q3 dq1 dq2 dq3] to [q1 dq1 q2 dq2 q3 dq3]
%     indqsource = 1:size(x0qdq, 2)/2;
%     inddqsource = ((size(x0qdq, 2)/2)+1):size(x0qdq, 2);
%     inqtarget = 1:2:size(x0qdq, 2);
%     indqtarget = 2:2:size(x0qdq, 2);
%     
%     x0(:, inqtarget) = x0qdq(:, indqsource);
%     x0(:, indqtarget) = x0qdq(:, inddqsource);
%     x(:, inqtarget) = xqdq(:, indqsource);
%     x(:, indqtarget) = xqdq(:, inddqsource);
    
function [H1, H2]= getRecoveryMatrix(iocInstance, prevH1, prevH2, x, u, dt)    
    % Get dynamic derivatives wrt to control
    % Get dynamic derivatives wrt to state
    [fx, fu, px, pu] = iocInstance.getDerivativesNewObservation(x, u);
    df_dx = fx';
    df_du = fu';
    dp_dx = px';
    dp_du = pu';
    
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


% % % function [H1, H2]= getRecoveryMatrix(iocInstance, prevH1, prevH2, x0qdq, u0, xqdq, u, dt)
% % %     x0 = x0qdq;
% % %     x = xqdq;
% % % 
% % %     % Get features derivatives
% % %     [~, ~, px, pu] = iocInstance.getDerivativesNewObservation(x, u0); 
% % %     dp_du = pu';
% % %     dp_dx = px';
% % %     
% % %     % Get dynamic derivatives wrt to control
% % %     [~, fu, ~, ~] = iocInstance.getDerivativesNewObservation(x0, u0); 
% % %     df_du = fu';
% % %     % Get dynamic derivatives wrt to state
% % %     % Derivatives x1
% % %     [fx, ~, ~, ~] = iocInstance.getDerivativesNewObservation(x, u);
% % %     df_dx = fx';
% % %     
% % %     if (isempty(prevH1) && isempty(prevH2))
% % %         H1=df_du*dp_dx+dp_du;
% % %         H2=df_du*df_dx;
% % %     else
% % %         H1 = [prevH1+prevH2*dp_dx;
% % %               df_du*dp_dx+dp_du];
% % %         H2 = [prevH2*df_dx;
% % %               df_du*df_dx];
% % %     end
% % % end
