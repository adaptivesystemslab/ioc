function plot_batch_featureSet(currFileEntry, filepathInstance, featureSet, ...
    ekf_markerMask, ekf_eventhandler, ekf_markerTally, modelInstance, algorithmParam)
    if ~isempty(featureSet)
        % plot figures and save
        savePath = globalConstants_filepaths.filePrefixJointsCombined(filepathInstance, currFileEntry, algorithmParam);
        featureSet.figureJointData(savePath, modelInstance, 'angle');
        
        savePath = globalConstants_filepaths.filePrefixVelocitiesCombined(filepathInstance, currFileEntry, algorithmParam);
        featureSet.figureJointData(savePath, modelInstance, 'velocity');

%        savePath = globalConstants_filepaths.filePrefixJointsIndividual(filepathInstance, currFileEntry, algorithmParam);
%        featureSet.figureJointData(savePath, modelInstance, 'angleInd');
        
%        savePath = globalConstants_filepaths.filePrefixMeasurementsIndividual(filepathInstance, currFileEntry, algorithmParam);
%        featureSet.figureMeasurements(savePath, ekf_markerMask, ekf_eventhandler, ekf_markerTally);

        switch algorithmParam.modalityType
            case 'imu'

            case 'mocap'
                savePath = globalConstants_filepaths.filePrefixMatchMatrix(filepathInstance, currFileEntry, algorithmParam);
                featureSet.figureMatchMatrix(savePath, modelInstance);
        end
    end
end