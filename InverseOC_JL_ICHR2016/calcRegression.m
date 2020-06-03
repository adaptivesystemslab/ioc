function [b, yCalc, Rsq] = calcRegression(x_in, y_in)
    X = [ones(length(x_in),1) x_in];
    b = X\y_in;
    yCalc = X*b;
    
    Rsq = 1 - sum((y_in - yCalc).^2)/sum((y_in - mean(y_in)).^2);
    
    if b(2) > 0
        % positive, do nothing
    else
        Rsq = -Rsq;
    end
end