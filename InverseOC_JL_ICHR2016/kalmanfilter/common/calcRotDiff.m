function dist = calcRotDiff(R1, R2)
    if size(R1, 2) == 4
        R1 = R1(1:3, 1:3);
        R2 = R2(1:3, 1:3);
    end
    
    dist_tmp = eye(3) - R1*R2';
    dist = 1 - (norm(dist_tmp, 'fro') / (2*sqrt(2))); % the divided factor is a normalizer so dist is now a percentage
end