function [cost_functions_used, corrValArray, corrValArrayOut, corrValId] = ...
    checkCorrelation(J_cost_all, J_const_all, costInd, functionMode, param)

loopFlag = 1;
cost_functions_used = 1:length(J_cost_all);

functionStr = strsplit(functionMode, '_');

while loopFlag
    J_cost_use_combined = [J_cost_all{:} J_const_all{:}];
    corrValArrayzz = normalizeMatrix(J_cost_use_combined);

    lenJvert = size(J_cost_use_combined, 1);
    lenJcost = length(J_cost_all);
    lenJconst = length(J_const_all);

    switch functionStr{2} 
        case '111'
             J_cost_use_combined_final = J_cost_use_combined;
             
        case '100'
            % if using 100, then report only the c/h results
            
            % keep only the q grad (vert) and q const (horz)
            indVert = 1:lenJvert/3;
            indHorz = [1:lenJcost (1:(lenJconst/3)) + lenJcost];

%             indVert = costInd; 
%             indHorz = [costInd costInd+lenJcost];
            J_cost_use_combined_final = J_cost_use_combined(indVert, indHorz);
            
%         case '010'
% %             indVert = (lenJvert/3+1):(2*lenJvert/3);
%             indVert = costInd + lenJvert/3;
%             J_cost_use_combined_final = J_cost_use_combined(indVert, :);
%             
%         case '001'
% %             indVert = (2*lenJvert/3 + 1):(3*lenJvert/3);
%             indVert = costInd + 2*lenJvert/3;
%             J_cost_use_combined_final = J_cost_use_combined(indVert, :);
    end
    
    switch functionStr{1}
        case 'corr'
            corrValArray = corrcoef(J_cost_use_combined_final);
    
        case 'dotp'
            corrValArray = dotProdCorr(J_cost_use_combined_final);
    end
    
    corrValArrayOut = triu(corrValArray, 1);
    corrValNonZero = corrValArrayOut ~= 0;
    [row, col] = find(corrValNonZero);
    rowcol = [row col];
    corrValArrayOut = corrValArrayOut(corrValNonZero);
    
    % also, remove the nans
    corrValNonNan = ~isnan(corrValArrayOut);
    corrValArrayOut = corrValArrayOut(corrValNonNan);
    rowcol = rowcol(corrValNonNan, :);
    corrValId = rowcol(:, 1) + 0.01*rowcol(:, 2);
    
    % check if any entries are above the thresold
    [maxVal, maxInd] = max(abs(corrValArrayOut));
    if ~isempty(corrValArrayOut) && maxVal > param.corrThreshold && param.removeHighlyCorrelatedJ
        fprintf('High correlation detected (%f > threshold %f), removing entry %u\n', maxVal, param.corrThreshold , maxInd);
        entryToRemove = max(rowcol(maxInd, :));
        J_cost_all{entryToRemove} = zeros(size(J_cost_all{1})); % keep it as zeros to keep the correlation matrix consistent
        cost_functions_used(entryToRemove) = 0;
    else
        loopFlag = 0;
    end
end
end

function corrValArray = normalizeMatrix(mtx)
% take the dot product of each col with each other
    corrValArray = zeros(size(mtx));
    for i1 = 1:size(mtx, 2)
        vec1 = mtx(:, i1) / norm(mtx(:, i1));
        corrValArray(:, i1) = vec1;
    end
end

function corrValArray = dotProdCorr(mtx)
% take the dot product of each col with each other
    corrValArray = zeros(size(mtx, 2), size(mtx, 2));
    for i1 = 1:size(mtx, 2)
        for i2 = 1:size(mtx, 2)
            vec1 = mtx(:, i1) / norm(mtx(:, i1));
            vec2 = mtx(:, i2) / norm(mtx(:, i2));

%             vec1 = mtx(:, i1) / max(abs(mtx(:, i1)));
%             vec2 = mtx(:, i2) / max(abs(mtx(:, i2)));
            
            corrValArray(i1, i2) = dot(vec1, vec2);
        end
    end
end