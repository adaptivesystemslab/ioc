function dx = calcDerivVert(x, dt)
%     dx = diff(x, 1, 1)/dt; 
%     dx = dx([1:end end], :);

    dx = [zeros(1, size(x, 2)); ...
        diff(x, 1, 1)/dt];