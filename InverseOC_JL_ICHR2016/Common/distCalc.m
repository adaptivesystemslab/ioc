function absDist = distCalc(ptA, ptB)
    % Given 2 points, determine the absolute distance between the two
    % handles 1D, 2D and 3D data
    
    arrayDim = size(ptA, 2);
    
    switch arrayDim
        case 1
            absDist = sqrt( ptA .^ 2 + ptB .^ 2 );
            
        case 2
            X1 = ptA(:,1);
            Y1 = ptA(:,2);
            
            X2 = ptB(:,1);
            Y2 = ptB(:,2);
            
            distX = X1 - X2;
            distY = Y1 - Y2;
            
            absDist = sqrt( distX .^ 2 + distY .^ 2 );
            
        case 3
            X1 = ptA(:,1);
            Y1 = ptA(:,2);
            Z1 = ptA(:,3);
            
            X2 = ptB(:,1);
            Y2 = ptB(:,2);
            Z2 = ptB(:,3);
            
            distX = X1 - X2;
            distY = Y1 - Y2;
            distZ = Z1 - Z2;
            
            absDist = sqrt( distX .^ 2 + distY .^ 2  + distZ .^ 2);
            
        otherwise
            absDist = 0;
    end
end