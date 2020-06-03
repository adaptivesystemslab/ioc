function pathToRawData = callDatasetBasePath(specStruct)
    % this function exist to serve a single modification point for all the
    % local installations for all data, as oppose to having to update the
    % pathing information in numerous different places. 

    basepath = 'C:\Documents\aslab\data\'; % if all your data is in one location, just update the basepath
    
    % update filepath to local installations
    switch lower(specStruct.datasetName) % move all to smallcaps
        case 'healthy1'
            pathToRawData = [basepath 'Lowerbody_healthy1_2011-11\'];

        case 'healthy2'
            pathToRawData = [basepath 'Lowerbody_healthy2_2013-07\']; % path to the raw data
            
        case 'stjoseph'
            pathToRawData = [basepath  'Lowerbody_StJoseph1_2013-02\'];
            
        case 'tri'
            pathToRawData = [basepath  'Lowerbody_TRI1_2012-10\'];
    end
    
end