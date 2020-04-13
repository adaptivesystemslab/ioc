function [weights_mean_all, weights_var_all, winCount_all] = cumWeights(t, progressVar, featureLabels)
    lenWeights = length(featureLabels);
    
    weights_mean_all = zeros(length(t), lenWeights);
    weights_var_all = zeros(length(t), lenWeights);
    winCount_all = zeros(length(t), 1);
    
    % the weight at each timestep is the sum of every window that overlaps with iterate
    for i = 1:length(t)
        weightAtI = [];
        for j = 1:length(progressVar)
            if isempty(progressVar(j).winInds)
                continue
            end
            
            if progressVar(j).winInds(1) <= i && progressVar(j).winInds(end) >= i
                if sum(isnan(progressVar(j).weights)) == 0
                    weightAtI = [weightAtI; progressVar(j).weights];
                end
            end
        end
        
        if ~isempty(weightAtI)
            weights_mean_all(i, :) = mean(weightAtI, 1);
            weights_var_all(i, :) = std(weightAtI, 0, 1);
            winCount_all(i, :) = size(weightAtI, 1);
        else
            weights_mean_all(i, :) = zeros(1, lenWeights);
            weights_var_all(i, :) = zeros(1, lenWeights);
            winCount_all(i, :) = 0;
        end
    end
end