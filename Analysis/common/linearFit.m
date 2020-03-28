function [b, Rsq2, X, yCalc2] = linearFit(x, y)
    X = [ones(length(x),1) x];
    b = X\y;
    yCalc2 = X*b;
    Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
end