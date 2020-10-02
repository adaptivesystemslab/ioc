function fct = setNormErrorFct()
%     rmseFct = @(x, y) sqrt( sum(sum((x - y).^2, 2)) / (size(x, 2)) );
%     rmseFct = @(x, y) sqrt( sum((normVector(x) - normVector(y)).^2) / (size(x, 1)) );
    fct = @(x, y) abs(normVector(x - y));