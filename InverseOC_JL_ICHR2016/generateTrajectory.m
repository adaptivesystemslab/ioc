function [q, dq, ddq] = generateTrajectory(angControl_y, angControl_x, t, angControl_dy, angControl_ddy, param)
    Nh = floor(length(t)/2);
    Nf = length(t);
    dofs = param.dof_count;
        
    % expect angControl points to be in arrays. if there is only one cell
    % array, then generate a second, identical, one
    if isempty(angControl_y)
        angControl_y = zeros(dofs, 2);
        angControl_dy = zeros(dofs, 2);
        angControl_ddy = zeros(dofs, 2);
        angControl_x(1) = 1;
        angControl_x(2) = Nf;
    elseif size(angControl_y, 2) == 1
        angControl_y(:, 2) = angControl_y(:, 1);
        angControl_dy(:, 2) = angControl_dy(:, 1);
        angControl_ddy(:, 2) = angControl_ddy(:, 1);
        angControl_x(1) = 1;
        angControl_x(2) = Nf;
    end
    
    % now evenly breakdown the time array to allow for even spacing between
    % all the control points
    indToUse_window = [];
%     win_shift = floor(length(t)/length(angControl_y));
    ind = 0;
    for ind_windowCount = 1:length(angControl_x) - 1
        ind = ind+1;
        
%         [closestVal1, closestInd1] = findClosestValue(angControl_x(1), t);
%         [closestVal2, closestInd2] = findClosestValue(angControl_x(2), t);
        
        indToUse_window(ind, 1) = angControl_x(ind_windowCount);
        indToUse_window(ind, 2) = angControl_x(ind_windowCount+1);
        
        if indToUse_window(ind, 2) > length(t)
            %         indToUse_window(ind, 2) = length(t_full);
            indToUse_window = indToUse_window(1:ind-1, :); % remove that latest entry
            break % we're at the end of the line
        end
    end
    
    if isempty(ind_windowCount)
        indToUse_window(1, 1) = 1; % forcable set the control point at the end of t
        indToUse_window(1, 2) = length(t); % forcable set the control point at the end of t
    else
        indToUse_window(ind, 2) = length(t); % forcable set the control point at the end of t
    end
    
    for j = 1:size(angControl_y, 2)-1
        % [X, dX, ddX, dddX] = Septic(q0, dq0, ddq0, dddq0, qf, dqf, ddqf, dddqf, t)
        indToUse = indToUse_window(j, 1):indToUse_window(j, 2);
        t_curr = t(indToUse);
        
        for k = 1:dofs
            sp_qd0 = quintic_poly(angControl_y(k, j),angControl_dy(k, j),angControl_ddy(k, j),...
                angControl_y(k, j+1),angControl_dy(k, j+1),angControl_ddy(k, j+1),t_curr);
            sp_qd1 = polyder(sp_qd0);
            sp_qd2 = polyder(sp_qd1);
            
            q(k, indToUse) =     polyval(sp_qd0, t_curr);
            dq(k, indToUse) =    polyval(sp_qd1, t_curr);
            ddq(k, indToUse) =   polyval(sp_qd2, t_curr);
        end
    end
end