function R = Rodrigues(w,dt)
    % w should be column vector of size 3

    % theta
    th = norm(w*dt);
    
%     wn = w;
    if (norm(w) == 0)
%         in case norm(w) zero and we're dividing by zero
        wn = w;
    else
        wn = w/norm(w);		% normarized vector
    end

    S = skewSym3(wn);
    R = eye(3) + S * sin(th) + S^2 * (1-cos(th));
end