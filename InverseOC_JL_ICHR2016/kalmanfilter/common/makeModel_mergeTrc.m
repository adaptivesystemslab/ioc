function modelStruct = makeModel_mergeTrc(modelInstance, trcData_all, pExternalParam, filepathModelXml)
    baseFrame = 'frame_6dof_root';
    baseMarker = 'HIP_BASE_JC';

    dataInstance_trc = rlDataInstance_trc(modelInstance);
    algorithmParam = structAlgorithmParam(pExternalParam);
    
    field = fieldnames(trcData_all);
    time = 1:length(trcData_all.(field{1}));
    dataInstance_trc.updateData(time, trcData_all);
    dataInstance_trc.dataProcessing(algorithmParam);
    
    modelInstance.model.base = baseFrame;
    dataInstance_trc.updateBaseMarker(baseMarker);
    
    [kinematicTransform, dynamicTransform, sensorTransform] = ...
        modelInstance.makeModel(filepathModelXml, dataInstance_trc, algorithmParam);
    
    % generate an 'average' marker frame that only has completed frames
    dataCurr = dataInstance_trc.measurement.getMesArray;
    dataToKeep = [];
    for i = 1:size(dataCurr, 1)
        missingInd = find(dataCurr(i, :) > algorithmParam.missingMarkerValue - 5);
        
        if isempty(missingInd)
            dataToKeep = [dataToKeep; dataCurr(i, :)];
        else
            lala = 1;
        end
    end
    
    modelStruct.kinematicTransform = kinTransOut;
    modelStruct.dynamicTransform = dynTransOut;
    modelStruct.sensorTransform = sensTransOut;
    modelStruct.sensorSecondaryTransform = sensSecTransOut;
    
    
    if 0 
        % load the model
        
    end
end
 