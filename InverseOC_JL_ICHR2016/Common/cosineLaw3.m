function gamma = cosineLaw3(A, B, C)
    % providing the three sides of a triangle, return the angle
    % COSINE LAW is c^2 = a^2 + b^2 - 2*a*b*cosd(gamma)
    % rearranged... gamma = arccos[ (a^2 + b^2 - c^2) / 2ab ]
    % -- Extended to 3D version (NOTE uncertain about accuracy)
    
    % --> Gamma is the angle where point 'C' is, or the angle opposite from
    % the length 'c'
    
    % for the purposes of KIN612L6 (left here for example purposes)
    % A = hand     -> a = upper arm
    % B = shoulder -> b = lower arm
    % C = elbow    -> c = shoulder to hand
    % gamma = elbow angle
    
    % -- DISTANCE --
%     a = distCalc(B, C);
%     b = distCalc(A, C);
%     c = distCalc(A, B);
    
    Xa = A(:, 1); Xb = B(:, 1); Xc = C(:, 1);
    Ya = A(:, 2); Yb = B(:, 2); Yc = C(:, 2);
    Za = A(:, 3); Zb = B(:, 3); Zc = C(:, 3);
    
    a = sqrt((Xc-Xb).^2 + (Yc-Yb).^2 + (Zc-Zb).^2);
    b = sqrt((Xc-Xa).^2 + (Yc-Ya).^2 + (Zc-Za).^2);
    c = sqrt((Xa-Xb).^2 + (Ya-Yb).^2 + (Za-Zb).^2);
    
    
    % -- ANGLE --
    precos = ( a .^ 2 + b .^ 2 - c .^ 2 ) ./ ( 2 * a .* b);
    gamma = acosd(precos);
end
