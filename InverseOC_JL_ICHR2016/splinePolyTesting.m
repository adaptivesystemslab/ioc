% load a function, plot pull out a window, and perform both spline and
% 5th order polynomial fit and report on its error
trajToLoadPath = '..\data\Subject01\Session1\SQUA_STD_NON1\Kinematics_meas_Sb1_Tr1.mat';
dynToLoadPath =  '..\data\Subject01\Session1\SQUA_STD_NON1\Dynamics_meas_Sb1_Tr1.mat';

param = setup_main(trajToLoadPath, dynToLoadPath, 'win');
param.h_vals =  [0];
param.h_array = param.h_vals*param.h; 

% % load the full trajectory
traj_load = load(trajToLoadPath);

t_full = param.t_full;
q_full = traj_load.q(2:4, :);
dq_full = traj_load.dq(2:4, :);
ddq_full = traj_load.ddq(2:4, :);

ind = 0;
for ind_windowCount = 1:param.win_shift:length(t_full)
    ind = ind+1;
    indToUse_window(ind, 1) = ind_windowCount;
    indToUse_window(ind, 2) = ind_windowCount + param.win_length - 1;
    
    if indToUse_window(ind, 2) > length(t_full)
%         indToUse_window(ind, 2) = length(t_full);
        indToUse_window = indToUse_window(1:ind-1, :); % remove that latest entry
        break % we're at the end of the line
    end
end
windowCount = size(indToUse_window, 1);

for ind_windowCount = 1:windowCount
    fprintf('(%u/%u) Loading window...\n', ind_windowCount, windowCount);
    
    % pull out the window of data
    q_opt = q_full(:, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2));
    dq_opt = dq_full(:, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2));
    ddq_opt = ddq_full(:, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2));
    t_opt = t_full(:, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2));
    
    % update the constraints based on windowed data
    const_y{1} = q_opt(:, 1  );
    const_y{2} = q_opt(:, 50 );
    const_y{3} = q_opt(:, 100);
    param.const_y = const_y;
    
    param.spline_endcond_x = [1 100 1 100];
    param.spline_endcond_y = [dq_opt(:, 1) dq_opt(:, 100) ddq_opt(:, 1) ddq_opt(:, 100)];
    
    % perform spline fit
    q_opt_knot = q_opt(:, param.intermed_ind);
    splineFit = calc_direct_spline(q_opt_knot, param);
    [q_spline, dq_spline1] = calc_features(splineFit, param);
            dq_spline2 = calcDeriv(q_spline, param.dt);
    rmse_spline(ind_windowCount) = sqrt(sum(sum((q_opt - q_spline) .^2, 2))/(size(q_opt, 1)));
                
    % perform poly fit
    [q_poly, dq_poly, polyParam] = calc_trajectory_polynomial(t_opt, q_opt);
    rmse_poly(ind_windowCount) = sqrt(sum(sum((q_opt - q_poly) .^2, 2))/(size(q_opt, 1)));
        
    % plot all data
    clf
    ax(1) = subplot(211);
    plot(t_full, q_full'); hold on
    title(['Spline: ' num2str(rmse_spline(ind_windowCount)) ', Poly: ' num2str(rmse_poly(ind_windowCount))]);   
    plot(t_opt, q_spline, 'o', 'DisplayName', 'Spline');
    plot(t_opt, q_poly,   'x', 'DisplayName', 'Poly');
    legend show
    
    ax(2) = subplot(212);
    plot(t_full, dq_full'); hold on
    plot(t_opt, dq_spline1, 'o', 'DisplayName', 'Spline1');
    plot(t_opt, dq_spline2, '+', 'DisplayName', 'Spline1');
    plot(t_opt, dq_poly,   'x', 'DisplayName', 'Poly');
    
    linkaxes(ax, 'x');
    
    xlim([t_opt(1)-0.2 t_opt(end)+0.2]);
end