function obj = arsLoader(varargin)
% arsLoader is a factory function that parses and interprets the header
% path that is passed in and generates the proper datahandle object
% accordingly, and load the data from file
% 
%  Usage:
%    obj = arsLoader(filepathHeader) will strip out the '.header'
%    extension from the filepathHeader and attach a '.data' to it for
%    filepathData
%      filepathHeader - header filepath
%
%    obj = arsLoader(filepathHeader, filepathData) performs the
%    same action as the first usage case, but the user can specify the 
%      filepathHeader - header filepath
%      filepathData - data filepath

    % parse the varargins
    if length(varargin) > 1
        filepathHeader = varargin{1};
        filepathData = varargin{2};
    elseif length(varargin) == 1
        % is it a .header path or a .data path?
        splitHeaderStr = strsplit(varargin{1},'.header');
        splitDataStr = strsplit(varargin{1},'.data');
        
        if length(splitHeaderStr) > 1
            % the variable passed in is a header
            filepathHeader = varargin{1};
            filepathData = [splitHeaderStr{1} '.data'];
        else
            filepathData = varargin{1};
            filepathHeader = [splitDataStr{1} '.header'];
        end
    else
        % unknown number of varargin
    end

    % load some computer specific settings. we currently don't have an easy
    % way to load XML using sharcnet, so bypass header loading for now if
    % using some flags
    switch computer
        case 'GLNXA64'
            % ignore the header files for SHARCNET, for now
            loadParam.loadHeader = 0;
            loadParam.loadData = 1;
            
        otherwise
            [idum,hostname] = system('hostname');
            
            if strcmpi(hostname, 'ROBOTLAB8')
                % only for JL's local machine
                loadParam.loadHeader = 1;
                loadParam.loadData = 1;
            else
                % load normal stuff for everyone else
                loadParam.loadHeader = 1;
                loadParam.loadData = 1;
            end
    end
    
    % load the header (XML file)
    if exist(filepathHeader, 'file') && loadParam.loadHeader
        xDoc = xmlread(filepathHeader);
        xRoot = xDoc.getDocumentElement;
        headerType = char(xRoot.getNodeName);
    elseif exist(filepathHeader, 'file') && ~loadParam.loadHeader
        % parse it straight up
        fId = fopen(filepathHeader);
        for i = 1:2
            tline = fgets(fId);
            headers = regexp(tline,' ','split');
            headerType = headers{1}(2:end);
            
            if ~strcmpi(headerType, '?xml')
                % not the string we're looking for. sometimes fgets seem to
                % ignore the first xml line. not sure why it doesn this
                break
            end
        end
        fclose(fId);
    elseif ~exist(filepathHeader, 'file') && loadParam.loadHeader && ...
            ~exist(filepathData, 'file') && loadParam.loadData
        % no files actually exist
        obj = [];
        return
        
    else
        % if there is no header, expect it to be EKF (for now)
        headerType = 'EKF';
    end
    
    % select the proper class based on header content
    switch headerType
        case 'EKF'
            % loading EKF data
            obj = ekfDataHandle(filepathHeader, filepathData, loadParam);
              
        case 'MotionCapture'
            % loading motion capture data
            obj = mocapDataHandle(filepathHeader, filepathData, loadParam);
            
        case 'Segmentation'
            % loading segmentation data
            obj = segmentationDataHandle(filepathHeader, filepathData, loadParam);
            
        case 'SensorPkg'
            % loading IMU data
            obj = imuDataHandle(filepathHeader, filepathData, loadParam);
            
        case 'Patient'
            % loading IMU data
            obj = patientDataHandle(filepathHeader, loadParam);
    end
    
    obj = obj.load(loadParam);
end