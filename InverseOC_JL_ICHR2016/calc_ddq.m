function [dq, ddq, dddq] = calc_ddq(q, param)
    dq = calcDeriv(q, param.dt_spline);
    ddq = calcDeriv(dq, param.dt_spline);
    dddq = calcDeriv(ddq, param.dt_spline);
end