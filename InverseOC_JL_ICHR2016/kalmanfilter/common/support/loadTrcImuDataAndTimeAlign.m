function [dataInstance_trc, dataInstance_imu] = loadTrcImuDataAndTimeAlign(currFileEntry, filepathTimeAlignment)
    % load mocap data
    modelInstance = rlModelInstance_ssa(0);
    dataInstance_trc = rlDataInstance_trc(modelInstance);
    dataInstance_trc.loadData(currFileEntry.filePathTrc);
    
    % load imu data
    dataInstance_imu = rlDataInstance_imu(modelInstance);
    dataInstance_imu.loadData(currFileEntry);
    
    % perform time alignment
    [trc, trcTime, imus, h] = timeAlign(dataInstance_trc.data, dataInstance_trc.time, dataInstance_imu.data);
    
    if ~isempty(h)
        saveas(h, filepathTimeAlignment, 'fig');
        saveas(h, filepathTimeAlignment, 'png');
        close(h);
    end
    
    % load the new time aligned data into dataInstance format
    dataInstance_trc.updateData(trcTime, trc);
    dataInstance_trc.dataProcessing();
    
    dataInstance_imu.updateData(trcTime, imus);
    
%     dataInstance_trc = [];
%     dataInstance_imu = [];
end