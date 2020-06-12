function [param, skip_count] = update_intermed_ind(param, knot_count)

if isfield(param, 'intermed_ind_set')
    param.intermed_ind = param.intermed_ind_set;
    skip_count = 0;
%     param.x_knot = param.x_knot_set;
else
    length_traj = size(param.t_spline, 2);
    skip_count = floor(length_traj/knot_count);
    intermed_ind = unique([1:skip_count:length_traj-1 length_traj param.const_x]); % this ensures there is a knot point at the end, and on constraints

    param.intermed_ind = intermed_ind;
end

param.x_knot = param.t_spline(param.intermed_ind); % update the x knots to ensure the dt is correct
