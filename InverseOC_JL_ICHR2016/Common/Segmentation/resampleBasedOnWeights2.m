function newDInd = resampleBasedOnWeights(D)
    % not sure how to do this. if it has a heavier weight,
    % include it more often, as par Seiffert2008
    
    % note that the returned newDInd is a function of the input D vector,
    % which may not be the full length of the original training data
    
    % pick 1/(n+1+wp) amounts of new points to keep, based on
    % the weights. change the wp as needed
    weightPenalty = 1;
        
    DUnique = flipud(unique(D)); %
    DInd = [];
    
    % generate gauss length
    gaussArray = gauss(0, length(DUnique), [1:length(DUnique)]');
    normalizedGaussArray = gaussArray/gaussArray(1);
    
    % expand out the normalizedGaussian
    normalizedGaussExpanded = zeros(size(D));
    for i = 1:length(DUnique)
        normalizedGaussExpanded(D == DUnique(i)) = normalizedGaussArray(i);
    end
    
    thresholdRand = rand(size(D));
    threshold = normalizedGaussExpanded > thresholdRand;
    newDInd = find(threshold);
    
%     
%     
%     % figure out how many entries are in each 'stage'
%     DCount = zeros(size(DUnique));
%     for ind_Dweights = 1:length(DUnique)
%         DCount(ind_Dweights) = length(find(D == DUnique(ind_Dweights)));
%     end
%     
%     % take the entry with the highest weights, and calculate resampling
%     % ratios from that one
%     countHeaviest = DCount(end);
%     DUse = zeros(size(DUnique));
%     for ind_Dweights = 1:length(DUnique)
%         denom = length(DUnique)-ind_Dweights+1+weightPenalty;
%         toUse = floor(countHeaviest/denom);
%         
%         if toUse > DCount(ind_Dweights)
%             % in case there is not enough entries in the lower weights
%             toUse = DCount(ind_Dweights);
%         end
%         
%         DUse(ind_Dweights) = toUse;
%         
%         newD = find(D == DUnique(ind_Dweights));
%         resampleInd = randperm(DUse(ind_Dweights)); % randomly pick some indices
%         newDToUse = newD(resampleInd(1:DUse(ind_Dweights)));
% 
% %         for ind_rep = 1:ind_Dweights
%             DInd = [DInd; newDToUse];
% %         end
%     end
% 
%     newDInd = sort(DInd);
end