function fct = setMavFct()
%     fct = @(x, y) sqrt( sum(sum((x - y).^2, 2)) / (size(x, 2)) );
%     fct = @(x, y) sqrt( sum((normVector(x) - normVector(y)).^2) / (size(x, 1)) );
    fct = @(x, y) ( sum((normVector(x - y))) / (size(x, 1)) );