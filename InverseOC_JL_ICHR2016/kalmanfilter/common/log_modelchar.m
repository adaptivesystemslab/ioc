function log_modelchar(modelInstance, filepath, currFileEntry, algorithmParam)
    if isempty(modelInstance)
        [fileID, errmsg] = fopen(filepath, 'a');
        fprintf(fileID,'%s,%s,%s',currFileEntry.subjectString, currFileEntry.exerciseName, '');
        fprintf(fileID,'%s','error');
        fprintf(fileID, '\n');
        fclose(fileID);
        return
    end

    if ~exist(filepath, 'file')
        %             if doesn't exist write header
        header = 'Subject,Exercise,ExerciseId,linkDef';
        
        for j = 1:length(modelInstance.kinematicTransform)
            header = [header ',normkin_' modelInstance.kinematicTransform(j).frameName ''];
        end
        header = [header ','];
%         for j = 1:length(modelInstance.dynamicTransform)
%             header = [header ',dyn_' modelInstance.dynamicTransform(j).name];
%         end
        for j = 1:length(modelInstance.sensorTransform)
            header = [header ',normsen_' modelInstance.sensorTransform(j).sensorName ''];
        end
        header = [header ','];
        for j = 1:length(modelInstance.sensorSecondaryTransform)
            header = [header ',normsensec_' modelInstance.sensorSecondaryTransform(j).sensorName ''];
        end
        header = [header ','];
        for j = 1:length(modelInstance.kinematicTransform)
            header = [header ',kin_' modelInstance.kinematicTransform(j).frameName '_x'];
            header = [header ',kin_' modelInstance.kinematicTransform(j).frameName '_y'];
            header = [header ',kin_' modelInstance.kinematicTransform(j).frameName '_z'];
        end
%         for j = 1:length(modelInstance.dynamicTransform)
%             header = [header ',dyn_' modelInstance.dynamicTransform(j).name];
%         end
        for j = 1:length(modelInstance.sensorTransform)
            header = [header ',sen_' modelInstance.sensorTransform(j).sensorName '_x'];
            header = [header ',sen_' modelInstance.sensorTransform(j).sensorName '_y'];
            header = [header ',sen_' modelInstance.sensorTransform(j).sensorName '_z'];
        end
        for j = 1:length(modelInstance.sensorSecondaryTransform)
            header = [header ',sensec_' modelInstance.sensorSecondaryTransform(j).sensorName '_x'];
            header = [header ',sensec_' modelInstance.sensorSecondaryTransform(j).sensorName '_y'];
            header = [header ',sensec_' modelInstance.sensorSecondaryTransform(j).sensorName '_z'];
        end        

        header = [header '\n'];
    else
        header = '';
    end
    
    % open/create log file
    [fileID, errmsg] = fopen(filepath, 'a');
    fprintf(fileID, header);
    
    % write data to log file
    if isnumeric(currFileEntry.fileId)
        fileId = num2str(currFileEntry.fileId);
    else
        fileId = currFileEntry.fileId;
    end
    fprintf(fileID,'%s,%s,%s,%s',currFileEntry.subjectString, currFileEntry.exerciseName, fileId, algorithmParam.linkDefinition);
    
    for j = 1:length(modelInstance.kinematicTransform)
        t = modelInstance.kinematicTransform(j).t;
        if isempty(t)
            t = nan(4, 4);
        end
        fprintf(fileID,',%.4f', norm(t(1:3, 4)'));
    end
    fprintf(fileID,',');
    for j = 1:length(modelInstance.sensorTransform)
        t = modelInstance.sensorTransform(j).t;
        if isempty(t)
            t = nan(4, 4);
        end
        fprintf(fileID,',%.4f', norm(t(1:3, 4)'));
    end
    fprintf(fileID,',');
    for j = 1:length(modelInstance.sensorSecondaryTransform)
        t = modelInstance.sensorSecondaryTransform(j).t;
        if isempty(t)
            t = nan(4, 4);
        end
        fprintf(fileID,',%.4f', norm(t(1:3, 4)'));
    end
    fprintf(fileID,',');
    for j = 1:length(modelInstance.kinematicTransform)
        t = modelInstance.kinematicTransform(j).t;
        if isempty(t)
            t = nan(4, 4);
        end
        fprintf(fileID,',%.4f', t(1:3, 4)');
    end
    for j = 1:length(modelInstance.sensorTransform)
        t = modelInstance.sensorTransform(j).t;
        if isempty(t)
            t = nan(4, 4);
        end
        fprintf(fileID,',%.4f', t(1:3, 4)');
    end
    for j = 1:length(modelInstance.sensorSecondaryTransform)
        t = modelInstance.sensorSecondaryTransform(j).t;
        if isempty(t)
            t = nan(4, 4);
        end
        fprintf(fileID,',%.4f', t(1:3, 4)');
    end

    fprintf(fileID, '\n');
    
    fclose(fileID);
end
