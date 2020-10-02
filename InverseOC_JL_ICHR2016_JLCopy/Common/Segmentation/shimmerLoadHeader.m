function [shimmerInfo] = shimmerLoadHeader(currExerPath, shimmerDevicesFiles)
    for i = 1:length(shimmerDevicesFiles)
        targetFile = fullfile(currExerPath, ['SensorData_' shimmerDevicesFiles{i} '.header']);
        doc = xmlread(targetFile);

        xRoot = doc.getDocumentElement;
        name = char(xRoot.getAttribute('name'));
        btId = char(xRoot.getAttribute('btId'));

        shimmerInfoTemp.name = name;
        shimmerInfoTemp.btId = btId;
        
        targetFile = fullfile(currExerPath, ['SensorData_' shimmerDevicesFiles{i} '.data']);
        shimmerData = parseCSV_shimmer(targetFile);

        shimmerInfoTemp.data = shimmerData;
        
        shimmerInfo{i} = shimmerInfoTemp;
    end
end