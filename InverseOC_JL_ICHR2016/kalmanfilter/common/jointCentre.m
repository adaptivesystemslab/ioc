function outmarkersArray = jointCentre(markers, distribution)
    if ~exist('distribution', 'var')
        per = 1/length(markers);
        for i = 1:length(markers)
            distribution{i} = [per per per];
        end
    end
    
    masking = ones(size(markers{1}));
    for i = 1:length(markers)
        normVec = normVector(markers{i});
        localMask = find(normVec == 0);
        if ~isempty(localMask)
            masking(localMask, :) = NaN;
        end
    end

    outmarkersArray = zeros(size(markers{1}));
    for i = 1:length(markers)
        distributionArray = repmat(distribution{i}, size(markers{i}, 1), 1);
        outmarkersArray = outmarkersArray + distributionArray.*masking.*markers{i};
    end
%     outmarkersArray = outmarkersArray / length(markers);
end