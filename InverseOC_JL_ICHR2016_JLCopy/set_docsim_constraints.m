function [param_docsim, set_docsim_check_struct] = set_docsim_constraints(param_docsim)
    % adding in velo/accel end conditions to ensure smooth curves
    switch param_docsim.doc_sim_sp_dq_endcond
        case 'constraint_dx'
            param_docsim.spline_endcond_dx = param_docsim.const_x;
            param_docsim.spline_endcond_dy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_dx));
%             param_docsim.spline_endcond_ddx = param_docsim.const_x;
%             param_docsim.spline_endcond_ddy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_ddx));
%             param_docsim.spline_endcond_dddx = param_docsim.const_x;
%             param_docsim.spline_endcond_dddy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_dddx));

        case 'two_dx'
            param_docsim.spline_endcond_dx = [1 param_docsim.spline_length];
            param_docsim.spline_endcond_dy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_dx));
%             param_docsim.spline_endcond_ddx = [1 param_docsim.spline_length];
%             param_docsim.spline_endcond_ddy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_ddx));
%             param_docsim.spline_endcond_dddx = [1 param_docsim.spline_length];
%             param_docsim.spline_endcond_dddy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_dddx));

        case 'two_dxddx'
            param_docsim.spline_endcond_dx = [1 param_docsim.spline_length];
            param_docsim.spline_endcond_dy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_dx));
            param_docsim.spline_endcond_ddx = [1 param_docsim.spline_length];
            param_docsim.spline_endcond_ddy = zeros(param_docsim.dof_count, length(param_docsim.spline_endcond_ddx));

        case 'none'
            param_docsim.spline_endcond_dx = [];
            param_docsim.spline_endcond_ddx = [];
            param_docsim.spline_endcond_dy = [];
            param_docsim.spline_endcond_ddy = [];
    end
    
    switch param_docsim.doc_sim_cf_dq_constraint
        case 'two_dx'
            param_docsim.const_dx = [1 param_docsim.spline_length];
            param_docsim.const_dy = zeros(param_docsim.dof_count, length(param_docsim.const_dx));
            param_docsim.const_ddx = [];
            param_docsim.const_ddy = [];

        case 'two_dxddx'
            param_docsim.const_dx = [1 param_docsim.spline_length];
            param_docsim.const_dy = zeros(param_docsim.dof_count, length(param_docsim.const_dx));
            param_docsim.const_ddx = [1 param_docsim.spline_length];
            param_docsim.const_ddy = zeros(param_docsim.dof_count, length(param_docsim.const_ddx));

        case 'constraint_dx'
            param_docsim.const_dx = param_docsim.const_x;
            param_docsim.const_dy = zeros(param_docsim.dof_count, length(param_docsim.const_dx));
            param_docsim.const_ddx = [];
            param_docsim.const_ddy = [];

        case 'constraint_dxddx'
            param_docsim.const_dx = param_docsim.const_x;
            param_docsim.const_dy = zeros(param_docsim.dof_count, length(param_docsim.const_dx));
            param_docsim.const_ddx = param_docsim.const_x;
            param_docsim.const_ddy = zeros(param_docsim.dof_count, length(param_docsim.const_ddx));

        case 'none'
            param_docsim.const_dx = [];
            param_docsim.const_dy = [];
            param_docsim.const_ddx = [];
            param_docsim.const_ddy = [];
    end
    
    % list the doc constraints
    lala = update_intermed_ind(param_docsim, param_docsim.doc_sim.knot_count);
    set_docsim_check_struct.intermed_ind = lala.intermed_ind;
    set_docsim_check_struct.const_x = param_docsim.const_x;
    set_docsim_check_struct.const_y = param_docsim.const_y;
    set_docsim_check_struct.const_dx = param_docsim.const_dx;
    set_docsim_check_struct.const_dy = param_docsim.const_dy;
    set_docsim_check_struct.const_ddy = param_docsim.const_ddy;
    set_docsim_check_struct.const_ddy = param_docsim.const_ddy;
    
    set_docsim_check_struct.spline_endcond_dx = param_docsim.spline_endcond_dx;
    set_docsim_check_struct.spline_endcond_dy = param_docsim.spline_endcond_dy;
    set_docsim_check_struct.spline_endcond_ddx = param_docsim.spline_endcond_ddx;
    set_docsim_check_struct.spline_endcond_ddy = param_docsim.spline_endcond_ddy;
    
    set_docsim_check_struct
end