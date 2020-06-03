function [tau, dtau, ddtau, cop, dcop, ddcop] = calc_dynamics_timeseries(q, dq, ddq, param)
    % calculate torque
    [tau, en] = SQUAT_7DOF_dyn3(q,dq,ddq,param.base0,param.ext_wrenches,param);

    F = en(1:3, :);
    M = en(4:6, :);
%         mappingJac = [param.COM1 param.COM2 param.COM3 param.COM4 param.COM5 param.COM6 param.COM7];
%         FP = FP2Base0([0 0 0], FP, mappingJac);
    %         grf = FP.N0(2:4, :);

    dtau = calcDeriv(tau, param.dt_spline);
    ddtau = calcDeriv(dtau, param.dt_spline);
    
    cop = CoP_calculation(en);
    dcop = calcDeriv(cop, param.dt_spline);
    ddcop = calcDeriv(dcop, param.dt_spline);
end

function FP = FP2Base0(base0, FP, mappingJac)
    F1=FP.F;
    N1=FP.M;

    r = mappingJac' - repmat(base0, 7, 1);

    for ii_fpcalc = 1:size(r, 1) % iterate through the joints
        r_inst = repmat(r(ii_fpcalc, :), size(F1, 2), 1);
        N0 = N1' + cross(r_inst, F1');
        N_out(ii_fpcalc, :) = N0(:, 3); % in our applications, only z rotation exists
    end

    FP.F0=F1;
    FP.N0=N_out;
end