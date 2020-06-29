function midArray = midpointCalc(ptA, ptB)
    % given two points, find the middle point

    arrayLen = size(ptA, 1);
    arrayDim = size(ptA, 2);
    
    midArray = zeros(size(ptA));

    for i = 1:arrayLen
        switch arrayDim
            case 1
                midArray(i) = average(ptA, ptB);

            case 2
                X1 = ptA(i,1);
                Y1 = ptA(i,2);

                X2 = ptB(i,1);
                Y2 = ptB(i,2);

                distX = average(X1, X2);
                distY = average(Y1, Y2);

                midArray(i, :) = [distX distY];

            case 3
                X1 = ptA(i,1);
                Y1 = ptA(i,2);
                Z1 = ptA(i,3);

                X2 = ptB(i,1);
                Y2 = ptB(i,2);
                Z2 = ptB(i,3);

                distX = average(X1, X2);
                distY = average(Y1, Y2);
                distZ = average(Z1, Z2);

                midArray(i, :) = [distX distY distZ];

            otherwise
                midArray = [];
                return
        end
    end
end

function val = average(ptA, ptB)
    val = mean([ptA ptB]);
end