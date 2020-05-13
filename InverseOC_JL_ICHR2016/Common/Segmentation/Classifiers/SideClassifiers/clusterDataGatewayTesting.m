function trainingDataParse = clusterDataGatewayTesting(trainingData, incsvmParam, dataToKeep)
    % note that for PP, an offset is applied if the whole array is being
    % loaded to line up with the stacking. but usually this function is
    % called when we're only looking at a window, so the offset will not be
    % triggered. apply the offset explicitly from the outside

    switch incsvmParam.dataFormatName
        case 'PC'
            if ~exist('dataToKeep', 'var')
                numberOfDataToKeep = floor(length(trainingData.testingLabel) * incsvmParam.percentageOriginalDataToKeepFirstReplacement);
                dataToKeep = sort(randperm(length(trainingData.testingLabel), numberOfDataToKeep));
                %                 dataToKeep = 1:length(trainingData.trainingLabel);
            elseif dataToKeep == 0
                dataToKeep = 1:length(trainingData.testingLabel);
            end

            trainingDataParse.input = trainingData.testingData(dataToKeep, :); % pre-PCA data
            trainingDataParse.output = trainingData.testingLabel(dataToKeep);
            trainingDataParse.outputGndTruth = trainingData.testingLabel(dataToKeep);

        case 'PP'
            if ~exist('dataToKeep', 'var')
                numberOfDataToKeep = floor(length(trainingData.label) * incsvmParam.percentageOriginalDataToKeepFirstReplacement);
                dataToKeep = sort(randperm(length(trainingData.label), numberOfDataToKeep));
                %                 dataToKeep = 1:length(trainingData.label);
            elseif dataToKeep == 0
%                 dataToKeep = 1+trainingData.settings.dimStrack:length(trainingData.label)-trainingData.settings.dimStrack;
                dataToKeep = dataToKeepTemp(trainingData.testingInd);
            end

            trainingDataParse.input = trainingData.data(dataToKeep, :); % the original data
            trainingDataParse.output = trainingData.label(dataToKeep);
            trainingDataParse.outputGndTruth = trainingData.label(dataToKeep);
    end