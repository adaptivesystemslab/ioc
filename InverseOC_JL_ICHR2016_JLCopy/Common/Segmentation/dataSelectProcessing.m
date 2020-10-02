function [data] = dataSelectProcessing(dataSelect, sigDofArrayUpdate)
    % select the data for training and testing purposes
    
%     dataSelect.settings.label_notSegment = 0;
%     dataSelect.settings.label_segment = 1;
    
    % if data has already been loaded (ie from a previous version 
    if strcmpi(dataSelect.settings.mode, 'training') && isfield(dataSelect, 'existingTrainingData') ...
            && ~isempty(dataSelect.existingTrainingData)
        data = dataSelect.existingTrainingData;
        
    else
        switch dataSelect.name
            case 'General'
                data = dataGeneral(dataSelect);
        end
    end
    
    if exist('sigDofArrayUpdate', 'var')
        data.updateSigDofArray(sigDofArrayUpdate);
    end
        
    data.partition;
end