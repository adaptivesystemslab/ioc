function [weights_mean_rank, weights_var_rank, winCount_rank] = cumWeights_rankFilt(t, progressVar, featureLabels, rankThres)
    lenWeights = length(featureLabels);
    
    weights_mean_rank = zeros(length(t), lenWeights);
    weights_var_rank = zeros(length(t), lenWeights);
    winCount_rank = zeros(length(t), 1);
    
    for i = 1:length(t)
        weightAtI = [];
        for j = 1:length(progressVar)
            if isempty(progressVar(j).winInds)
                continue
            end
            
            if progressVar(j).winInds(1) <= i && progressVar(j).winInds(end) >= i && progressVar(j).rankTraj(end) >= rankThres
                weightAtI = [weightAtI; progressVar(j).weights];
            end
        end
        
        if size(weightAtI, 1) > 0
            weights_mean_rank(i, :) = mean(weightAtI, 1);
            weights_var_rank(i, :) = std(weightAtI, 1);
            winCount_rank(i, :) = size(weightAtI, 1);
        else
            weights_mean_rank(i, :) = zeros(1, lenWeights);
            winCount_rank(i, :) = 0;
        end
    end
end