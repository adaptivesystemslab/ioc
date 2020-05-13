function [trainingData] = pullDataForCircleMagic(currTime, currData, windowSize, coeffCount)
    
%     windowSize = 15;
%     coeffCount = 6;
    [a b] = butter(5,[0.1/128*2 10/128*2]);
    
    currTime = currTime - currTime(1);    
    currData = filtfilt(a,b,currData);
    dof = size(currData, 2);

    if size(currData, 1) == 1
        trainingInd = 1;
        windowSize = 1;
    else
        trainingInd = windowSize+1:size(currData, 1)-windowSize-1;
    end
    
    options = statset('TolFun',5e-3,'MaxIter',30);

    for ind_segPt = 1:length(trainingInd)
        % pull out the entries for the 'segment' examples
        currInd = trainingInd(ind_segPt);
        startInd = currInd-windowSize;
        endInd = currInd+windowSize;

        % copy out the array, and reshape it so it'll stack into the
        % training matrix
        currWindowedTime = currTime(startInd:endInd) - currTime(startInd);
        currWindowedData = currData(startInd:endInd, :);
        currWindowedReshape = reshape(currWindowedData, [1 dof*(windowSize*2+1)]);

        coeffs = nlinfit(currWindowedTime,currWindowedReshape,...
            @modelfun,[ones(1,coeffCount) zeros(1,coeffCount)],options); % 6x2 coeff for the regression
        
        trainingData(ind_segPt, :) = coeffs;
    end
end