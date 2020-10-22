function dx = calcDeriv(x, dt)
    % row vector
    dx = diff(x, 1, 2)/dt; 
    dx = dx(:, [1:end end]);
