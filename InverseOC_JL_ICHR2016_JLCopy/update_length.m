function param = update_length(param, win_length, spline_length)
% abstract out the length changing components of the code
% - win_length
% - spline_length

    if nargin == 2
        spline_length = win_length;
    end

    param.win_length = win_length;
    param.spline_length = spline_length;
    
    param.dt_spline = param.dt_full * (win_length/spline_length);
    
    param.t_full = param.ti_full:param.dt_full:param.tf_full; % full traj
    param.t_spline = (param.ti_full:(param.spline_length-1))*param.dt_spline; % just the window, what we will be calc'ing over
    param.NbSample = length(param.t_spline); 
end

