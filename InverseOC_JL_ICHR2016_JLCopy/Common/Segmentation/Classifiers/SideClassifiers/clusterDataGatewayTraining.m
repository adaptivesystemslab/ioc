function trainingDataParse = clusterDataGatewayTraining(trainingData, incsvmParam, dataToKeep)
    % note that for PP, an offset is applied if the whole array is being
    % loaded to line up with the stackin
    
    % if it's not empty, load the testing data
%     switch mode
%         case 'testing'
%             switch incsvmParam.dataFormatName
%                 case 'PC'
%                     if ~exist('dataToKeep', 'var')
%                         numberOfDataToKeep = floor(length(trainingData.testingLabel) * incsvmParam.percentageOriginalDataToKeepFirstReplacement);
%                         dataToKeep = sort(randperm(length(trainingData.testingLabel), numberOfDataToKeep));
%                         %                 dataToKeep = 1:length(trainingData.trainingLabel);
%                     elseif dataToKeep == 0
%                         dataToKeep = 1:length(trainingData.testingLabel);
%                     end
% 
%                     trainingDataParse.input = trainingData.tData(dataToKeep, :); % pre-PCA data
%                     trainingDataParse.output = trainingData.testingLabel(dataToKeep);
%                     trainingDataParse.outputGndTruth = trainingData.testingLabel(dataToKeep);
% 
%                 case 'PP'
%                     if ~exist('dataToKeep', 'var')
%                         numberOfDataToKeep = floor(length(trainingData.label) * incsvmParam.percentageOriginalDataToKeepFirstReplacement);
%                         dataToKeep = sort(randperm(length(trainingData.label), numberOfDataToKeep));
%                         %                 dataToKeep = 1:length(trainingData.label);
%                     elseif dataToKeep == 0
%                         dataToKeep = 1+trainingData.settings.dimStrack:length(trainingData.label)-trainingData.settings.dimStrack;
% 
%                     end
% 
%                     trainingDataParse.input = trainingData.data(dataToKeep, :); % the original data
%                     trainingDataParse.output = trainingData.label(dataToKeep);
%                     trainingDataParse.outputGndTruth = trainingData.label(dataToKeep);
%             end
%             
%         case 'training'
            switch incsvmParam.dataFormatName
                case 'PC'
                    if ~exist('dataToKeep', 'var')
                        numberOfDataToKeep = floor(length(trainingData.trainingLabel) * incsvmParam.percentageOriginalDataToKeepFirstReplacement);
                        dataToKeep = sort(randperm(length(trainingData.trainingLabel), numberOfDataToKeep));
                        %                 dataToKeep = 1:length(trainingData.trainingLabel);
                    elseif dataToKeep == 0
                        dataToKeep = 1:length(trainingData.trainingLabel);
                    end

                    trainingDataParse.input = trainingData.trainingData(dataToKeep, :); % pre-PCA data
                    trainingDataParse.output = trainingData.trainingLabel(dataToKeep);
                    trainingDataParse.outputGndTruth = trainingData.trainingLabel(dataToKeep);
                    trainingDataParse.sigdofArray = trainingData.trainingSigdofcellarray(dataToKeep);

                case 'PP'
                    if ~exist('dataToKeep', 'var')
                        numberOfDataToKeep = floor(length(trainingData.trainingInd) * incsvmParam.percentageOriginalDataToKeepFirstReplacement);
                        dataToKeepVal = sort(randperm(length(trainingData.trainingInd), numberOfDataToKeep));
                        dataToKeep = trainingData.trainingInd(dataToKeepVal);
                        %                 dataToKeep = 1:length(trainingData.label);
                    elseif dataToKeep == 0
%                         dataToKeep = (1+trainingData.settings.dimStrack):(length(trainingData.label)-trainingData.settings.dimStrack);
                        dataToKeep = trainingData.trainingInd;
                    end
                    


                    trainingDataParse.input = trainingData.data(dataToKeep, :); % the original data
                    trainingDataParse.output = trainingData.label(dataToKeep);
                    trainingDataParse.outputGndTruth = trainingData.label(dataToKeep);
                    trainingDataParse.sigdofArray = trainingData.sigdofcellarray(dataToKeep);
            end
%     end