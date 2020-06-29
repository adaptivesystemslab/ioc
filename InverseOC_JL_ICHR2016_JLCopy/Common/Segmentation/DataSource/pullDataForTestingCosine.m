function [trainingData, trainingLabel] = pullDataForTraining(obj, windowSize)
   
    segmentPt = 1:length(obj.label);
    
    trainingData = [];
    trainingLabel = [];
    
    for ind_segPt = 1:length(segmentPt)
        % pull out the entries for the 'segment' examples
        currInd = segmentPt(ind_segPt);
        startInd = currInd-windowSize;
        endInd = currInd+windowSize;
        
        if startInd < 1 % skip the edge cases (easiest to deal with 
            continue
%             startInd = 1;
        end
        
        if endInd > length(obj.label) % skip the edge cases (easiest to deal with 
            continue
%             endInd = length(obj.label);
        end

        % copy out the array, and reshape it so it'll stack into the
        % training matrix
        currWindowedData = obj.data(startInd:endInd);
        currWindowedReshape = reshape(currWindowedData, [1 windowSize*2+1]);

        trainingData = [trainingData; currWindowedReshape];
        trainingLabel = [trainingLabel; obj.label(currInd)];
    end
    
%     trainingLabel = ones(size(trainingData, 1), 1)*labelType;
end