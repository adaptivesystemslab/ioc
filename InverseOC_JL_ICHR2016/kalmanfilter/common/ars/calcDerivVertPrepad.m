function dx = calcDerivVertPrepad(x, dt)
%     dx = diff(x, 1, 1)/dt; 
%     dx = dx([1:end end], :);

    diffVal = diff(x, 1, 1)/dt;
    
    dx = [diffVal(1, :); ...
        diffVal];