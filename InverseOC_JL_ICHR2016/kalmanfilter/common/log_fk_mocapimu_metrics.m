function measRmseId1 = calc_fk_mocapimu_metrics(featureSet_fk, featureSet_orig, filepathFig, filepathImuFkMocapIndividualLog, currFileEntry,externalParam_orig, ...
    ekf_markerMask, ekf_eventhandler, ekf_markerTally)

    [measRmseId1, measRmseId2, sensorsUsed, sensorsTotal] = ...
        plot_measurements(featureSet_orig, featureSet_fk, featureSet_orig.measurement, featureSet_fk.measurement, 'orig', 'fk', filepathFig, ...
        ekf_markerMask, ekf_eventhandler, ekf_markerTally);
 
    meanMeasRmseId1 = mean(measRmseId1);
    meanMeasRmseId2 = mean(measRmseId2);
    sensorsUsedPercentage = sum(sensorsUsed)/sum(sensorsTotal);
    
    % now write some stuff to file
    if ~exist(filepathImuFkMocapIndividualLog, 'file')
        %             if doesn't exist write header
        header = 'Subject,Exercise,linkDef';
        header = [header ',meanRmseId1,meanRmseId2,sensorPercentage,'];
        
        for j = 1:length(featureSet_orig.measurement_labels)
            header = [header ',RMSE1_' featureSet_orig.measurement_labels{j}];
        end
        for j = 1:length(featureSet_orig.measurement_labels)
            header = [header ',RMSE2_' featureSet_orig.measurement_labels{j}];
        end

        header = [header '\n'];
    else
        header = '';
    end

    % open/create log file
    [fileID, errmsg] = fopen(filepathImuFkMocapIndividualLog, 'a');
    fprintf(fileID, header);
    
    % write data to log file
    fprintf(fileID,'%s,%s,%s',currFileEntry.subjectString, currFileEntry.exerciseName, externalParam_orig.linkDefinition );
    fprintf(fileID,',%f,%f,%f', meanMeasRmseId1, meanMeasRmseId2, sensorsUsedPercentage);
    
    for j = 1:length(featureSet_orig.measurement_labels)
        fprintf(fileID,',%.4f', measRmseId1(j));
    end
    for j = 1:length(featureSet_orig.measurement_labels)
        fprintf(fileID,',%.4f', measRmseId1(j));
    end

    fprintf(fileID, '\n');
    
    fclose(fileID);
end
