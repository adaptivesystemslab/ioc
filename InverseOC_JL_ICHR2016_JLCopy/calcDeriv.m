function dx = calcDeriv(x, dt)
    dx = diff(x, 1, 2)/dt; 
    dx = dx(:, [1:end end]);
