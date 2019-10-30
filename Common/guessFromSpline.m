function guess = guessFromSpline(model, initialState, finalState)
    [q, dq, tau, time] = getSpline(model, initialState, finalState);
           
    guess = struct('q', q, 'dq', dq, 'tau', tau, 'time', time);
end

