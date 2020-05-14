function  writeHeader_segment(targetFile, segmentParam, deviceInfo)
% make header data for segmentation. expects a few variables
%   targetFile - the full path the where the header should be made
%   segmentParam
%     .segMethod - (string) “cropping”, “manual”, “automatic”
%     .date - (datenum) date when the segmentation happened 
%     .segFacilitator - (string) any additional labels to how it is done would be here
%     .sourceData - (string) the data that was used for segmentation
%     .badSegments - (string) a string array of bad segment indices
%     .comments - (string) any additional comments left by the segmenter 
%   deviceInfo (optional) 
%     .deviceId - (cell of strings) btId or other identifier for the devices used
%     .deviceLocation - (cell of strings) the location to where the device is mounted
%     .initGyroVal - (optional) (cell of double) initial gyroscope value for DC offset
    
    writeGlobalHeader(targetFile);
    templateVariables = setTemplateVariables;
    
    fId = fopen(targetFile, 'w');
        
    fprintf(fId, '<Segmentation type="%s">\n', segmentParam.segMethod);
    fprintf(fId, '  <Facilitator date="%s" details="%s" sourceData="%s"/>\n', ...
        datestr(segmentParam.date, templateVariables.datestr), ...
        segmentParam.segFacilitator, segmentParam.sourceData);
    
    if exist('deviceInfo', 'var')
        fprintf(fId, '  <Devices>\n');
        
        for i = 1:length(deviceInfo.deviceId)
            fprintf(fId, '    <Device id="%s" location="%s">\n', ...
                deviceInfo.deviceId{i}, deviceInfo.deviceLocation{i});
            
            if isfield(deviceInfo, 'initGyroVal')
                fprintf(fId, '      <GyroOffsetRaw>%f,%f,%f</GyroOffsetRaw>\n', deviceInfo.initGyroVal{i});
            end
            
            fprintf(fId, '    </Device>\n');
        end
        
        fprintf(fId, '  </Devices>\n');
    end

    fprintf(fId, '  <comments badSegments="%s" comments="%s"/>\n', ...
        segmentParam.badSegments, segmentParam.comments);    
    fprintf(fId, '</Segmentation>\n');
    
    fclose(fId);
end