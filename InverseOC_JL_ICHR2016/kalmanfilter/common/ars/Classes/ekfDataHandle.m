classdef ekfDataHandle < ADataHandle
%  ekfDataHandle is a model to handle data parsing from EKF data files. 
%  Pass in the header and the data file path to the EKF data file and it 
%  will load both the header and the data structure 
%
%  Usage:
%    ekfData = ekfDataHandle(filepathHeader, filepathData) 
%      filepathHeader - header filepath
%      filepathData - data filepath
%
%  Example:
%    Loading data from file
%      basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy2_2013-07\Subject2\Session1\HAAO_STD\EKF\2014_02_13\';
%      headerPath = [basePath 'ekf.header'];
%      dataPath = [basePath 'ekf.data'];
%      ekfData = ekfDataHandle(headerPath, dataPath);
%      ekfData = ekfData.load; 
    
properties (SetAccess = public)
%     properties % the following are inheritated from ADataHandle
%         % handle information
%         filepathHeader = '';
%         filepathData = '';
%     end

    % header 
    modelPath = [];
    model;
    date = []; % date that the EKF was ran
    processNoise = [];
    observationNoise = [];

    % data 
    data = []; % data as is from the data file

    % data (reparsed from the data array)
    time = []; 
    dt = [];
    Q = []; % if the motion contains joint angle data, it is copied here
    dQ = [];
    ddQ = [];

    % data from the file structure that is currently not being stored:
    % FrameSeqNumCalibrated
    % TimeStampUncalibrated,TimeStampCalibrated,SystemMsTimeStampCalibrated,SystemTimeStampCalibrated,
end % properties

methods
    function obj = ekfDataHandle(varargin)
        % constructor function. parse the incoming variables
        obj.filepathHeader = varargin{1};
        obj.filepathData = varargin{2};
    end

    function obj = load(obj, varargin)
        % load the data from file based on header and data properties
        if nargin > 0
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
    
    function obj = importFromModel(obj, pNoise, oNoise, state, time, dt)
        obj.date = now; % date that the EKF was ran
        obj.processNoise = pNoise;
        obj.observationNoise = oNoise;
        
        obj.data = state;
        obj.time = time;
        obj.dt = ones(size(time))*dt;
    end

    function write(obj)
        % write the stored header and data to file
        obj.writeHeader;
        obj.writeData;
    end

    function h = plot(obj,varargin)
        h = figure;

        % zero offset the whole thing
        initTime = obj.time(1);
        time = obj.time - initTime; % a local instance of time, since we don't want to change the overall time array

%         h1 = subplot(2, 1, 1);
        plot(time, rad2deg(obj.Q));
        ylabel('Joint angle [deg]');

        if ~isempty(obj.model)
            jointStr = {obj.model.joints.name};
            
            for i = 1:length(jointStr)
                jointStr{i} = jointStr{i}(5:end);
            end
            legend(jointStr);
        end
        
%         h2 = subplot(2, 1, 2);
%         plot(time, rad2deg(obj.dQ));
%         ylabel('Joint velocity [deg/s]');

        xlabel('Time [s]');

%         linkaxes([h1, h2], 'x');
    end
    
    function cropData(obj, timeStart, timeEnd)
        [timeStartVal, timeStartInd] = findClosestValue(timeStart, obj.time);
        [timeEndVal, timeEndInd] = findClosestValue(timeEnd, obj.time);
        
        % apply the cropping
        obj.time = obj.time(timeStartInd:timeEndInd);
        obj.Q = obj.Q(timeStartInd:timeEndInd, :);
        obj.dQ = obj.dQ(timeStartInd:timeEndInd, :);
        obj.ddQ = obj.ddQ(timeStartInd:timeEndInd, :);
    end
    
    function loadModel(obj, folderpath)
        currModelPath = [folderpath '\' obj.modelPath];
        fprintf('Load model %s\n', currModelPath);
        obj.model = rlCModel(currModelPath);
    end
end % methods

methods (Access = private)
    function obj = loadHeader(obj, loadParam)

        % todo if the header doesn't exist, fill in some default values
        if ~exist(obj.filepathHeader, 'file')
            obj.date = now;
            return
        end
        
        % parse the XML file into a structure
        xDoc = parseXML(obj.filepathHeader);

        % pull the proper data from the base layer
        for ind_attributes = 1:length(xDoc.Attributes)
            xmlName = xDoc.Attributes(ind_attributes).Name;
            xmlValue = xDoc.Attributes(ind_attributes).Value;

            switch xmlName
                case 'date'
                    obj.date = datenum(xmlValue);
            end
        end

        for ind_attributes = 1:length(xDoc.Children)
            xmlName = xDoc.Children(ind_attributes).Name;
            xmlData = xDoc.Children(ind_attributes).Data;
            xmlChild = xDoc.Children(ind_attributes).Children;

            switch xmlName
                case '#text' 
                    if sum(isspace(xmlData)) > 0
                        % parsed a blank. do nothing. not sure why this
                        % happens
                    end

                case 'ProcessNoise'
                    obj.processNoise = eval(['[' xmlChild.Data '];']); % the child data is parsed as a string, so this converts it back to an array
                    
                    % process noise should be square, but may be written as
                    % a flat array
                    if size(obj.processNoise, 1) ~= size(obj.processNoise, 2)
                        obj.processNoise = diag(obj.processNoise);
                    end
                    
                case 'ObservationNoise'
                    obj.observationNoise = eval(['[' xmlChild.Data '];']);
                    
                    % obs noise should be square, but may be written as
                    % a flat array
                    if size(obj.observationNoise, 1) ~= size(obj.observationNoise, 2)
                        obj.observationNoise = diag(obj.observationNoise);
                    end
                    
                case 'model'
                    obj.modelPath = xmlChild(3).Children.Data;
            end
        end
    end
   
    function obj = loadData(obj, loadParam)
       % load the data via parseCSV
        ekfData = parseCSV(obj.filepathData);

        % check the time vector
        if isfield(ekfData, 'EKFTimeStamp')
            ekfData.SystemSTimeStamp = ekfData.EKFTimeStamp;
        else
            day2minMultiplier = 24*60*60;
            %     assessTimeMS = datenum('1970-1-1') + ekfData.SystemMSTimeStamp(1)/(1000*day2minMultiplier);
            assessTimeS = datenum('1970-1-1') + ekfData.SystemMSTimeStamp(1)/(day2minMultiplier);
            
            if assessTimeS > now
                % if this statement is true, then it suggests that the timestamp is
                % in milliseconds and should be scaled down by 1000
                ekfData.SystemSTimeStamp = ekfData.SystemMSTimeStamp/1000;
            else
                ekfData.SystemSTimeStamp = ekfData.SystemMSTimeStamp;
            end
        end
     
        fileType = 'MATLAB';
        if isfield(ekfData, 'DT')
            ekfData.Dt = ekfData.DT;
        elseif isfield(ekfData, 'DeltaTime')
            ekfData.Dt = ekfData.DeltaTime;
            fileType = 'ARS';
        end
        
        % now consolidate the joint angles
        if 0
            % if there some header info we can use, use those
        
        elseif isfield(ekfData, 'dQ1')
            % legacy
            qArray = [ekfData.Q1 ekfData.Q2 ekfData.Q3 ekfData.Q4 ekfData.Q5];
            dqArray =  [ekfData.dQ1 ekfData.dQ2 ekfData.dQ3 ekfData.dQ4 ekfData.dQ5];
            ddqArray =  [ekfData.ddQ1 ekfData.ddQ2 ekfData.ddQ3 ekfData.ddQ4 ekfData.ddQ5];
            
        elseif isfield(ekfData, 'QD1')
            % switch it over to the other type ('dQ' vs 'QD') - legacy
            % support
            ekfData.dQ1 = ekfData.QD1;
            ekfData.dQ2 = ekfData.QD2;
            ekfData.dQ3 = ekfData.QD3;
            ekfData.dQ4 = ekfData.QD4;
            ekfData.dQ5 = ekfData.QD5;

            ekfData.ddQ1 = ekfData.QDD1;
            ekfData.ddQ2 = ekfData.QDD2;
            ekfData.ddQ3 = ekfData.QDD3;
            ekfData.ddQ4 = ekfData.QDD4;
            ekfData.ddQ5 = ekfData.QDD5;
            
            fields = {'QD1','QD2','QD3','QD4','QD5','QDD1','QDD2','QDD3','QDD4','QDD5'};
            ekfData = rmfield(ekfData,fields);
            
            qArray = [ekfData.Q1 ekfData.Q2 ekfData.Q3 ekfData.Q4 ekfData.Q5];
            dqArray =  [ekfData.dQ1 ekfData.dQ2 ekfData.dQ3 ekfData.dQ4 ekfData.dQ5];
            ddqArray =  [ekfData.ddQ1 ekfData.ddQ2 ekfData.ddQ3 ekfData.ddQ4 ekfData.ddQ5];
            
        elseif isfield(ekfData, 'S1')
            % figure out how many states are available in total
            ekfFieldNames = fieldnames(ekfData);
            totalStateCount = -1; % state count starts at 0
            for ind_ekfFieldNames = 1:length(ekfFieldNames)
                currEkfFieldName = ekfFieldNames{ind_ekfFieldNames};
                if strcmpi(currEkfFieldName(1), 'S') && all(ismember(currEkfFieldName(2), '0123456789')) 
                    % attempt to parse the rest as a number
                    totalStateCount = totalStateCount + 1;
                end
            end
            
            % assume that the states divide evenly into 3: q, dq, ddq
            qStateCount = floor((totalStateCount+1)/3);
            
            qArray = [];
            dqArray = [];
            ddqArray = [];

            for ind_stateCount = (1:qStateCount)-1
                qState = ['S' num2str(ind_stateCount)];
                dqState = ['S' num2str(ind_stateCount+qStateCount)];
                ddqState = ['S' num2str(ind_stateCount+qStateCount*2)];
                
                % assign to individual dims
                eval(['ekfData.Q' num2str(ind_stateCount) ' = ekfData.' qState ';']);
                eval(['ekfData.dQ' num2str(ind_stateCount) ' = ekfData.' dqState ';']);
                eval(['ekfData.ddQ' num2str(ind_stateCount) ' = ekfData.' ddqState ';']);
                
                % assign to large matrix
                eval(['qArray = [qArray ekfData.Q' num2str(ind_stateCount) '];']);
                eval(['dqArray = [dqArray ekfData.dQ' num2str(ind_stateCount) '];']);
                eval(['ddqArray = [ddqArray ekfData.ddQ' num2str(ind_stateCount) '];']);
                
                % remove it from the original
                fields = {qState, dqState, ddqState};
                ekfData = rmfield(ekfData,fields);
            end
        end

        % saving the proper loaded data to the obj
        obj.time = ekfData.SystemSTimeStamp;
        obj.dt = ekfData.Dt;
        obj.Q = qArray;
        obj.dQ = dqArray;
        obj.ddQ = ddqArray;
    end

    function writeHeader(obj)
        %Write XML header
        docNode = com.mathworks.xml.XMLUtils.createDocument('EKF');
        docRootNode = docNode.getDocumentElement;
        docRootNode.setAttribute('date',datestr(now));
        
        %Process Noise
        el = docNode.createElement('ProcessNoise');
        str = sprintf('%f, ',diag(obj.processNoise));
        el.appendChild(docNode.createTextNode(str(1:end-2)));
        docRootNode.appendChild(el);
        
        %Observation Noise
        el = docNode.createElement('ObservationNoise');
        str = sprintf('%f, ',diag(obj.observationNoise));
        el.appendChild(docNode.createTextNode(str(1:end-2)));
        docRootNode.appendChild(el);
        
        xmlwrite(obj.filepathHeader,docNode);
    end

    function writeData(obj)
        
        writeDataStruct.SystemMSTimeStamp = obj.time;
        writeDataStruct.Dt = obj.dt;
        
        % parse all the states
        for i = 1:size(obj.data, 2)
            eval(['writeDataStruct.S' num2str(i) ' = obj.data(:, ' num2str(i) ');']);
        end
        
        writeCSV(obj.filepathData, writeDataStruct);
    end
    

end % methods (private)
end % classdef