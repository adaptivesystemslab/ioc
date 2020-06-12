classdef segmentationDataHandle < ADataHandle
%  segmentationDataHandle is a model to handle data parsing from segment
%  data files. Pass in the header and data file system to the seg file and 
%  it will load both  the header and the data structure. 
%
%  Usage:
%    ekfData = segmentationDataHandle(filepathHeader, filepathData) 
%      filepathHeader - header filepath
%      filepathData - data filepath
%
%  Example:
%    Loading data from file
%      basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy2_2013-07\Subject2\Session1\HAAO_STD\Segmentation_manual\';
%      headerPath = [basePath 'SegmentData_Manual_Manual.header'];
%      dataPath = [basePath 'SegmentData_Manual_Manual.data'];
%      ekfData = ekfDataHandle(headerPath, dataPath);
%      ekfData = ekfData.load; 

properties (SetAccess = public)
%     properties % the following are inheritated from ADataHandle
%         % handle information
%         filepathHeader = '';
%         filepathData = '';
%     end

    % header 
    segmentType = [];
    date = []; % date that the data collection occured
    details = []; % details regarding the segmentation process
    sourceData = []; % how were the segments created
    devices = []; % sampling rate of the IMU
    badSegments = []; % the type of sensor used
    comments = [];

    % data 
    data = []; % data as is from the data file

    % data (reparsed from the data array)
    use = [];
    use_train = [];
    use_test = [];
    segmentCount = [];
    timeStart = [];
    timeEnd = [];
    
    primitiveId = [];
    segmentId = [];    
    segmentName = [];
    extraFields = [];
    % data from the file structure that is currently not being stored:
    % FrameSeqNumCalibrated
    % TimeStampUncalibrated,TimeStampCalibrated,SystemMsTimeStampCalibrated,SystemTimeStampCalibrated,
end % properties

methods
    function obj = segmentationDataHandle(varargin)
        % constructor function. parse the incoming variables
        obj.filepathHeader = varargin{1};
        obj.filepathData = varargin{2};
    end
    
    function obj = updateFilePath(obj, varargin)
        % constructor function. parse the incoming variables
        obj.filepathHeader = varargin{1};
        obj.filepathData = varargin{2};
    end

%     function obj = load(obj)
%         % load the data from file based on header and data properties
%         obj = obj.loadHeader; % TODO - why is this not presistent (ie why does not obj.loadHeader work?)
%         obj = obj.loadData;
%     end
    
    function obj = load(obj, varargin)
        % load the data from file based on header and data properties
        if ~isempty(varargin)
            loadParam = varargin{1};
            
            if ~isfield(loadParam, 'loadHeader')
                loadParam.loadHeader = 1;
                loadParam.loadData = 1;
            end
        else
            loadParam.loadHeader = 1;
            loadParam.loadData = 1;
        end
        
        if loadParam.loadHeader
            obj = obj.loadHeader; % TODO - why is this not presistent (ie why does not obj.loadHeader work?)
        end
        
        if loadParam.loadData
            obj = obj.loadData;
        end
    end
    
    function obj = import(obj, varargin)
    
        p = inputParser;
        p.KeepUnmatched = true;
        
       addOptional(p,'segmentType','');
       addOptional(p,'date',now);
       addOptional(p,'details','');
       addOptional(p,'sourceData','');
       addOptional(p,'devices',[]);
       addOptional(p,'badSegments','');
       addOptional(p,'comments','');
       
       addOptional(p,'use',[]);
       addOptional(p,'use_train',[]);
       addOptional(p,'use_test',[]);
       addOptional(p,'timeStart',[]); % not really option, but double arrays fail on isnumeric in the parser
       addOptional(p,'timeEnd',[]);
       addOptional(p,'primitiveId',[]);
       addOptional(p,'segmentId',[]);
       addOptional(p,'segmentName',[]);
       addOptional(p,'extraFields',[]);
       
       parse(p,varargin{:});
       
       % insert into the headers
       obj = obj.importHeader(p.Results);
       obj = obj.importData(p.Results);
    end
    
    function write(obj)
        % write the stored header and data to file
        obj.writeHeader;
        obj.writeData;
    end

    function h = plot(obj, varargin)
        h = [];
        fprintf('Segmentation plots are non-sensical without the joint angle data. Please call the plotting function from the execise level. \n');
    end
    
    function struct = convertToStruct(obj)
        struct.use = obj.use;
        struct.segmentCount = obj.segmentCount;
        struct.timeStart = obj.timeStart;
        struct.timeEnd = obj.timeEnd;
    end
    
    function segInfo = extractSubsegmentation(obj, segId)
        % produce a new object that only has the segId
        primIndex = find(obj.primitiveId == segId);
        
        segInfo = retainSegments(obj, primIndex);
    end
    
    function segInfo = retainSegments(obj, indRetain)
        % produce a new object that only has the segId
        dataSegmentSubset.use = obj.use(indRetain);
        dataSegmentSubset.timeStart = obj.timeStart(indRetain);
        dataSegmentSubset.timeEnd = obj.timeEnd(indRetain);
        dataSegmentSubset.primitiveId = obj.primitiveId(indRetain);
        dataSegmentSubset.segmentId = obj.segmentId(indRetain);
        dataSegmentSubset.segmentName = obj.segmentName(indRetain);
        
        if ~isempty(obj.use_train)
            dataSegmentSubset.use_train = obj.use_train(indRetain);
        end
        
        if ~isempty(obj.use_train)
            dataSegmentSubset.use_test = obj.use_test(indRetain);
        end
        
        segInfo = segmentationDataHandle([], []);
        segInfo.import(dataSegmentSubset);
    end
    
    function segInfo = dropUseZeroSegments(obj)
        keepInds = find(obj.use == 1);
        segInfo = retainSegments(obj, keepInds);
    end
    
    function use_applyTrain(obj)
        obj.use = obj.use_train;
    end
    
    function use_applyTest(obj)
        obj.use = obj.use_test;
    end
    
    function obj = makeCopy(obj)
        
    end
end % methods

methods (Access = private)
    function obj = loadHeader(obj)

        % parse the XML file into a structure
        xDoc = parseXML(obj.filepathHeader);

        % pull the proper data from the base layer
        for ind_attributes = 1:length(xDoc.Attributes)
            xmlName = xDoc.Attributes(ind_attributes).Name;
            xmlValue = xDoc.Attributes(ind_attributes).Value;

            switch xmlName
                case 'type'
                    obj.segmentType = xmlValue;
            end
        end

        for ind_attributes = 1:length(xDoc.Children)
            xmlName = xDoc.Children(ind_attributes).Name;
            xmlData = xDoc.Children(ind_attributes).Data;
            xmlAttributes = xDoc.Children(ind_attributes).Attributes;
            xmlChild = xDoc.Children(ind_attributes).Children;
            
            switch xmlName
                case '#text'
                    if sum(isspace(xmlData)) > 0
                        % parsed a blank. do nothing. not sure why this
                        % happens
                    end
                    
                case 'Facilitator'
                    % read calibration data
                    for ind_childAttributes = 1:length(xmlAttributes) 
                        xmlChildName = xmlAttributes(ind_childAttributes).Name;
                        xmlChildValue = xmlAttributes(ind_childAttributes).Value;
                        
                        switch xmlChildName
                            case 'date'
                                obj.date = datenum(xmlChildValue, 'yyyy_mm_dd');
                            case 'details'
                                obj.details = xmlChildValue;
                            case 'sourceData'
                                obj.sourceData = xmlChildValue;
                        end
                    end
                    
                case 'Devices'
                    for ind_calibration = 1:length(xmlChild)
                        obj.devices{ind_calibration}.sensorId = xmlChild(ind_calibration).Attributes(1).Value;
                        obj.devices{ind_calibration}.name = xmlChild(ind_calibration).Attributes(2).Value;
                        obj.devices{ind_calibration}.gyroOffset = [0,0,0];
                        if ~isempty(xmlChild(ind_calibration).Children)
                           if ~isempty(xmlChild(ind_calibration).Children.Children)
                               if ~isempty(xmlChild(ind_calibration).Children.Children.Data)
                                  obj.devices{ind_calibration}.gyroOffset = str2num(xmlChild(ind_calibration).Children.Children.Data);
                               end
                           end
                        end
                    end
                    
                case 'comments'
                    for ind_childAttributes = 1:length(xmlAttributes)
                        xmlChildName = xmlAttributes(ind_childAttributes).Name;
                        xmlChildValue = xmlAttributes(ind_childAttributes).Value;
                        
                        switch xmlChildName
                            case 'badSegments'
                                obj.badSegments = xmlChildValue;
                            case 'comments'
                                obj.comments = xmlChildValue;
                        end
                    end
            end
        end
    end

    function obj = loadData(obj)
        % load the data via parseCSV. should never skip last entry
        segData = parseCSV(obj.filepathData, 0);
        
        % legacy renaming
        if isfield(segData, 'TimeStart')
            segData.timeStart = segData.TimeStart;
            segData.timeEnd = segData.TimeEnd;
            segData.use = segData.Use;
        end
        
        % check the time vector
        day2minMultiplier = 24*60*60;
        %     assessTimeMS = datenum('1970-1-1') + ekfData.SystemMSTimeStamp(1)/(1000*day2minMultiplier);
        assessTimeS = datenum('1970-1-1') + segData.timeStart(1)/(day2minMultiplier);
        
        if assessTimeS > now
            % if this statement is true, then it suggests that the timestamp is
            % in milliseconds (vlad style) and should be scaled down by 1000
            segData.timeStart =  segData.timeStart/1000;
            segData.timeEnd =    segData.timeEnd/1000;
        else
%             segData.timeStart =  segData.timeStart;
%             segData.timeEnd =    segData.timeEnd;
        end
    
        % saving the proper loaded data to the obj
        obj.use = segData.use;
        obj.segmentCount = 1:length(segData.use);
        obj.timeStart = segData.timeStart;
        obj.timeEnd = segData.timeEnd;
    end
    
    function obj = importHeader(obj,p)
        % load header from another source
        obj.segmentType = p.segmentType;
        obj.date = p.date; % date that the data collection occured
        obj.details = p.details; % details regarding the segmentation process
        obj.sourceData = p.sourceData; % how were the segments created
        obj.devices = p.devices; % sampling rate of the IMU
        obj.badSegments = p.badSegments; % the type of sensor used
        obj.comments = p.comments;
    end
    
    function obj = importData(obj,p)
        % load data from another source
        if isempty(p.use)
            p.use = ones(length(p.timeStart), 1);
        end
        
        obj.use = p.use;
        obj.timeStart = p.timeStart;
        obj.timeEnd = p.timeEnd;
        
        obj.primitiveId = p.primitiveId;
        obj.segmentId = p.segmentId;
        obj.segmentName = p.segmentName;
        obj.extraFields = p.extraFields;
        obj.use_train = p.use_train;
        obj.use_test = p.use_test;
        
        if isempty(obj.primitiveId)
            obj.primitiveId = zeros(size(obj.timeStart));
        end
        
        if isempty(obj.segmentId)
            obj.segmentId = zeros(size(obj.timeStart));
        end
        
        if isempty(obj.segmentName)
            obj.segmentName = cell(size(obj.timeStart));
        end
    end
        
    function writeHeader(obj)
        % make header data for segmentation
        writeGlobalHeader(obj.filepathHeader);
        templateVariables = setTemplateVariables;
        
        fId = fopen(obj.filepathHeader, 'a');
        
        fprintf(fId, '<Segmentation type="%s">\n', obj.segmentType);
        fprintf(fId, '  <Facilitator date="%s" details="%s" sourceData="%s"/>\n', ...
            datestr(obj.date, templateVariables.datestr), ...
            obj.details, ...
            obj.sourceData);
        
        fprintf(fId, '  <Devices>\n');
        for i = 1:length(obj.devices)
            fprintf(fId, '    <Device id="%s" location="%s">\n', obj.devices{i}.sensorId, obj.devices{i}.sensorName);
            if isfield(obj.devices, 'initGyroVal')
                fprintf(fId, '      <GyroOffsetRaw>%f,%f,%f</GyroOffsetRaw>\n', obj.devices{i}.initGyroVal);
            end
            fprintf(fId, '    </Device>\n');
        end
        fprintf(fId, '  </Devices>\n');
        
        fprintf(fId, '  <comments badSegments="%s" comments="%s"/>\n', ...
            obj.badSegments, obj.comments);
        fprintf(fId, '</Segmentation>\n');
        
        % write to the output
        fclose(fId);
    end

    function writeData(obj)
        % write data to file
        segmentDataToWrite.segmentCount = (1:length(obj.timeStart))';
        segmentDataToWrite.use = reshape(obj.use,[],1);
        segmentDataToWrite.timeStart = obj.timeStart;
        segmentDataToWrite.timeEnd = obj.timeEnd;
        writeCSV(obj.filepathData, segmentDataToWrite); % writeCSV expects all data array to be vertical
    end
end % methods (private)
end % classdef