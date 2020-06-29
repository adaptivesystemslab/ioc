function [trainingData, trainingLabel, trainingInd, trainingName] = pullDataForTrainingPatients(obj, sampleRules, labelType, windowSize, permissibleIndStruct, returnLimit)
    
    coeffCount = 6;
    options = statset('TolFun',5e-3,'MaxIter',30);

    % pull out the px as specified
    dataTemp = [];
    labelTemp = [];
    nameTemp = {};
%     jointVectorInd = obj.subjectData{1}.jointVectorInd;
    
    for ind_subject = obj.ind_partition
        currData = obj.subjectData{ind_subject};
        
        if isempty(sampleRules.px) || sum(currData.subjectNumber == sampleRules.px)
            dataTemp = [dataTemp; currData.data];
            labelTemp = [labelTemp; currData.label];
            nameTemp = [nameTemp; currData.name];
        end
    end
    
    % then select all the labelType points
    if ~isempty(labelType)
        trainingInd = find(labelTemp == labelType);
        
        % if specified, remove the first and last cluster if the labels
%         if sampleRules.cropFirstLastPoints
%             diffTrainingArray = [1; diff(trainingInd)];
%             gapTrainingArrayMid = find(diffTrainingArray > 1);
%             
%             % figure out the points to remove
%             pointsToRemoveStart = 1:gapTrainingArrayMid(1);
%             pointsToRemoveEnd = trainingInd(gapTrainingArrayMid(end)):trainingInd(end);
%             pointsToRemove = [pointsToRemoveStart pointsToRemoveEnd];
%             
%             trainingInd = setdiff(trainingInd, pointsToRemove)';
%         end
        
        % if specified, downsample the points away from the other motion type
        if sampleRules.tightSampling
%             trainingInd = codeSampleTightenThreshold(trainingInd, labelTemp);
            trainingInd = codeSampleTightenGaussian(trainingInd, labelTemp);
        end
    else
        trainingInd = 1:length(labelTemp);
    end
    
    dof = size(dataTemp, 2);
    trainingData = [];
    trainingLabel = [];
    
    % remove the points that are above or below the window size so error
    % checking is not needed later on
    lessInd = find(trainingInd-windowSize > 0);
    moreInd = find(trainingInd+windowSize <= length(labelTemp));
    intersectInd1 = trainingInd(intersect(lessInd, moreInd));
    
    % apply the velocity thresholds from permissibleIndStruct.veloThreshold
    intersectInd2 = intersect(permissibleIndStruct.veloThreshold, intersectInd1); 
    
    % apply the downsampling calculated from permissibleIndStruct.downSample
    intersectInd3 = intersect(permissibleIndStruct.downSample, intersectInd2);
    
    if sampleRules.cropFirstLastPoints
        intersectInd4 = intersect(permissibleIndStruct.SEP, intersectInd3);
    else
        intersectInd4 = intersectInd3;
    end
    
    % apply the p_1 limits from the data source
    intersectInd5 = intersect(permissibleIndStruct.segInclude, intersectInd4);
    
    intersectInd = intersectInd5;
    
    % if there is more points here than the 'return limit', truncate
    % trainingInd a bit
    if length(intersectInd) > returnLimit
        % randomly select some points
        selectInd = randperm(length(intersectInd), returnLimit);
        trainingInd = sort(intersectInd(selectInd));
    else
        trainingInd = intersectInd;
    end
      
    targetSize = [length(trainingInd) dof*(windowSize*2+1)]; % dof*(windowSize*2+1)+dof/2*coeffCount*2
%     targetSize = [length(trainingInd) 12]; % dof/2*coeffCount*2
    trainingData = zeros(targetSize);
    trainingLabel = zeros(targetSize(1), 1);
    trainingName = cell(size(trainingLabel));
    
    for ind_segPt = 1:length(trainingInd)
        % pull out the entries for the 'segment' examples
        currInd = trainingInd(ind_segPt);
        startInd = currInd-windowSize;
        endInd = currInd+windowSize;
        
%         if startInd < 1 % skip the edge cases (easiest to deal with 
%             continue
% %             startInd = 1;
%         end
%         
%         if endInd > length(labelTemp) % skip the edge cases (easiest to deal with 
%             continue
% %             endInd = length(obj.label);
%         end

        % copy out the array, and reshape it so it'll stack into the
        % training matrix     

        
        currWindowedCum = zeros(1, dof*windowSize*2+1);
        currWindowedRegression = zeros(1, targetSize(2));
        currTime = 0:1/128:(windowSize*2)/128;
        for ind_jointVector = 1:size(dataTemp, 2)
            cumStartInd = (ind_jointVector-1)*(windowSize*2+1)+1;
            cumEndInd = (ind_jointVector-1)*(windowSize*2+1)+windowSize*2+1;
            
            currWindowedData = dataTemp(startInd:endInd, ind_jointVector);
            currWindowedReshape = reshape(currWindowedData, [1 (windowSize*2+1)]);
            currWindowedCum(cumStartInd:cumEndInd) = currWindowedReshape;
            
%             % run regression fit on the first 5 (assuming it's q)
% %             if ind_jointVector > 0 && ind_jointVector < 6
%             if ind_jointVector == 4
%                 regStartInd = (ind_jointVector-1)*(coeffCount*2)+1;
%                 regEndInd = (ind_jointVector-1)*(coeffCount*2)+coeffCount*2;
%             
%                 coeffs = nlinfit(currTime,currWindowedReshape',...
%                     @modelfun,[ones(1,coeffCount) zeros(1,coeffCount)],options); % 6x2 coeff for the regression
%                 currWindowedRegression(:) = coeffs(:);
%             end
        end
        
%         currWindowedCum = [currWindowedRegression];

        trainingData(ind_segPt, :) = currWindowedCum;
        trainingLabel(ind_segPt) = labelTemp(currInd);
        
        if ~isempty(nameTemp)
            trainingName(ind_segPt) = nameTemp(currInd);
        end
        
%         if length(trainingLabel) >= returnLimit
%             return
%         end
    end
    
%     trainingLabel = ones(size(trainingData, 1), 1)*labelType;
end