function [A, Bn] = inverseCross(A, B, C, initGuess)
    % given A x B = C, this returns A
    % usage: A = inverseCross(A, B, C, accel, w, initGuess)
    % in context:
    % A = alpha
    % B = r
    % C = accel term in RNE
    
    % CONSTANTS
    upperBound = 0.04;
    lowerBound = 0.01;
    
    % works for motions in one direction
%     Bn = initGuess;
    
    % but we want flexible code, no?
   if (norm(A) < 10^-20)
%         make a good initial guess
        fprintf(['InverseCross: Calling normalizer \n']);
        Bn = initGuess;
    else
        AnormA = A/norm(A);

        dotP = dot(B, AnormA);
        Bp = (AnormA)*dotP;

        Bn = B - Bp;
%         Bn = bnCheck(B - Bp, upperBound, lowerBound);
    end

    projectedDotProduct = dot(A, Bn);
    if (abs(projectedDotProduct) > 10^-25) 
        % this should be zero
%         fprintf(['InverseCross: Projection is non-zero: ', num2str(projectedDotProduct), '\n']);
    else
%         fprintf(['InverseCross: Projection is good:     ', num2str(projectedDotProduct), '\n']);
    end

%     S = skewSym3(w);
%     C = accel - S*S*Bn; % dana's method
%     C = accel - S*S*B; % given way

% if frame == 2 && (norm(A) > 10^-25)
%     A'
%     [B'; Bp'; Bn']
%     C'
%     projectedDotProduct
% %     applepie = 5;
% end
    
    A = cross(Bn, C) / norm(Bn) ^ 2; % inverse cross
end

function BnNew = bnCheck(Bn, upper, lower)
    BnNew = zeros(size(Bn));
    
    for x = 1:length(Bn)
        BnExamine = abs(Bn(x));  % abs the elements
        
        if BnExamine > lower && BnExamine < upper
            % value is larger than the lower bound, but smaller than the
            % upper bound (so between the brackets), it is likely a
            % numerical error and should be supressed to zero
            BnNew(x) = 0;
            fprintf(['InverseCross: Zeroing numberical error:     ', num2str(BnExamine), '\n']);
        else
            BnNew(x) = BnExamine;
        end
    end
end