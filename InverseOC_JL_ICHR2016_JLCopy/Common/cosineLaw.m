function gamma = cosineLaw(A, B, C)
    % providing the three sides of a triangle, return the angle
    % COSINE LAW is c^2 = a^2 + b^2 - 2*a*b*cosd(gamma)
    % rearranged... gamma = arccos[ (a^2 + b^2 - c^2) / 2ab ]
    
    % --> Gamma is the angle where point 'C' is, or the angle opposite from
    % the length 'c'
    
    % for the purposes of KIN612L6 (left here for example purposes)
    % A = hand     -> a = upper arm
    % B = shoulder -> b = lower arm
    % C = elbow    -> c = shoulder to hand
    % gamma = elbow angle
    
    % -- DISTANCE --
    a = distCalc(B, C);
    b = distCalc(A, C);
    c = distCalc(A, B);
    
    % -- ANGLE --
    precos = ( a .^ 2 + b .^ 2 - c .^ 2 ) ./ ( 2 * a .* b);
    gamma = acosd(precos);
end

