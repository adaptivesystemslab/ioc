function spline_fit = calc_inverse_spline(y_knot_orig, param)
% generate the spline trajectory from the estimated path

len_h = length(param.h_array);
len_knot = length(param.intermed_ind);
len_dofs = param.dof_count;

spline_fit(len_dofs, len_knot, len_h).y_knot = []; % preallocation
spline_fit(len_dofs, len_knot, len_h).sp_qd0 = []; 
spline_fit(len_dofs, len_knot, len_h).sp_qd1 = []; 
spline_fit(len_dofs, len_knot, len_h).sp_qd2 = []; 
spline_fit(len_dofs, len_knot, len_h).sp_qd3 = []; 
for ii_h = 1:len_h % number of joints
    for ii_knot = 1:len_knot
        for ii_dofs = 1:len_dofs
            y_knot = y_knot_orig;

            % modify a knot point based on h, over each knot point. we want
            % to perturb each state (individual knot/dof combination) at a
            % time for the numerical gradient calculation
            y_knot(ii_dofs, ii_knot) = y_knot(ii_dofs, ii_knot) + param.h_array(ii_h);  % +, 0, or - h to each knot and calc the spline       
            splineFit = calc_direct_spline(y_knot, param);

            % dim: dof, knots, h
            spline_fit(ii_dofs, ii_knot, ii_h) = splineFit;
        end
    end
end
end
