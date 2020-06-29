function [q, dq, ddq, dddq, t_range] = calc_q_syms(splineFit, param)
    switch splineFit.type                   
        case 'piecewise_cds'
            for ind_n = 1:length(splineFit.sp_qd0) % note that pcds are constructed differently than the other spline types
                param_local = param;
%                 param_local.t_spline = splineFit.sp_qd0{ind_n}.t_spline;
                param_local.t_spline = splineFit.t_spline;
                
                for ind_m = 1:length(splineFit.sp_qd0)
                    q{ind_n}{ind_m} =     polyval_syms(splineFit.sp_qd0{ind_n}.sp_qd0{ind_m}, param_local.t_spline);
                    dq{ind_n}{ind_m} =    polyval_syms(splineFit.sp_qd0{ind_n}.sp_qd1{ind_m}, param_local.t_spline);
                    ddq{ind_n}{ind_m} =   polyval_syms(splineFit.sp_qd0{ind_n}.sp_qd2{ind_m}, param_local.t_spline);
                    dddq{ind_n}{ind_m} =  polyval_syms(splineFit.sp_qd0{ind_n}.sp_qd3{ind_m}, param_local.t_spline);
                end
                
                switch ind_n
                    case 1
                        ind_start = splineFit.sp_qd0{ind_n}.x_knot(1);
                        ind_end = splineFit.sp_qd0{ind_n}.x_knot(end);
                    otherwise
                        ind_start = splineFit.sp_qd0{ind_n}.x_knot(1)+1; % this keeps it causal
                        ind_end = splineFit.sp_qd0{ind_n}.x_knot(end);
                end
                
                 t_range(ind_start:ind_end) = ind_n;
            end
    end
end

function q_syms = polyval_syms(p,x)
    syms t q_syms
    for k = 1:length(p)
        q_syms(k) = p(k)*t^(k-1);
    end
    
    q_syms = sum(q_syms);
end