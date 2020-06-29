function trcOut = loadDataFromTrc_expressiveioc(pathToRawData, algorithmParam)

    trc = parseTRC(pathToRawData);        
    field_names = fieldnames(trc.data);
    
    if exist('algorithmParam', 'var') && ~isempty(algorithmParam.ekfRun)
        indsToKeep = algorithmParam.ekfRun;
    else
        indsToKeep = 1:length(trc.data.Time);
    end
    
    % Add a small value to RAC marker, since it cannot be zero
%     trc.data.RCAJ = trc.data.RCAJ + 1e-6;
    
    trcOut.time = trc.data.Time(indsToKeep);
    trcOut.dt = 1/trc.DataRate;
    
    % remove shoulder offset from every frame
%     offsetval = trc.data.RCAJ / 1000;
    
    % apply a world rotation to the mocap data
    offsetrot = rotx(90);
    %offsetrot = eye(3);
    
    for i=1:numel(field_names)
        sensorName = field_names{i};
        
        if(size(trc.data.(field_names{i}),2) == 3) % keeps only entries that are nX3, ie marker data    
            %trcOut.data.(sensorName) = (trc.data.(sensorName)')' / 1000; % data to [mm]
            temp = (trc.data.(sensorName)')' / 1000; % data to [mm]
%             trcOut.data.(sensorName) = trcOut.data.(sensorName) - offsetval;
            trcOut.data.(sensorName) = (offsetrot*temp(indsToKeep, :)')';
        end
    end
end