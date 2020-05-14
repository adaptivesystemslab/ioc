function d = matrixNorm(A, B, type)

d = 0;
switch type
    case 1
        for i = 1:size(A, 1)
            for j = 1:size(A, 2)
                d = d + abs(A(i, j) - B(i, j));
            end
        end
        
    case 2
        for i = 1:size(A, 1)
            for j = 1:size(A, 2)
                d = d + abs(A(i, j) - B(i, j))^2;
            end
        end
        d = sqrt(d);
        
    case Inf
        amax = max(max(A));
        bmax = max(max(B));
        d = abs(amax - bmax);
        
    case 'cos'
%         normA = sqrt(sum(A .^ 2, 2));
%         normB = sqrt(sum(B .^ 2, 1));
%         dtemp = bsxfun(@rdivide, bsxfun(@rdivide, A * B, normA), normB);
%         d = sum(sum(dtemp));

        cosdist = 0;
        for i = 1:size(A, 2)
            avec = A(:, i);
            bvec = B(:, i);
            cosvec = dot(avec, bvec) / (norm(avec)*norm(bvec));
            cosdist = cosdist + cosvec;
        end
        d = cosdist;
        
    case 'fro'
        dtemp = (A-B) * (A-B)';
        d = sqrt(trace(dtemp));
        
    otherwise
        
        
end