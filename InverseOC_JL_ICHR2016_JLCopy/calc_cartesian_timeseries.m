function [x, dx, ddx, dddx] = calc_cartesian_timeseries(q, param)

baseLink = [1;0;0];

x = zeros((7+2)*3, length(param.t));
for ii = 1:length(param.t)
    [Tr, x(:, ii)] = forward_kin(baseLink,q(:,ii),param);
end

% plotX(x, param)

% modify so it only keeps the 3 links, and discard base link and ankle link 
% pos since it's not moving. so ind 7 to ind 15
dx   = calcDeriv(x, param.dt);
ddx  = calcDeriv(dx, param.dt);
dddx = calcDeriv(ddx, param.dt);
end