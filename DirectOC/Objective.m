function f = Objective(t, x, u, weights, iocObject)

    w = weights; w = w/sum(w);
        
    % Compute current time increase and update value in iocInstance
    dt = t(2) - t(1);
    iocObject.dt = dt;
    
    features = iocObject.calcFeatures(x', u');
    f = sum(w.*features,2)';
    
end
