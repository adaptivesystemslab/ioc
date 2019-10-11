function guess = guessFromSpline(model, initialState, finalState)
    [q, dq, tau, time] = getSpline(model, initialState, finalState);
           
    guess = struct('jointAngles', q, 'angularVelocities', dq,...
               'control', tau, 'time', time);
end

