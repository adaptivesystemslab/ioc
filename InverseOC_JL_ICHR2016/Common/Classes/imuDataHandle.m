classdef imuDataHandle < ADataHandle
%  imuDatahandle is a model to handle data parsing from IMU data files. 
%  Pass in the header and data file path to the IMU file and it will load
%  both the header and the data structure 
%
%  Usage:
%    ekfData = imuDatahandle(filepathHeader, filepathData) 
%      filepathHeader - header filepath
%      filepathData - data filepath
%
%  Example:
%    Loading data from file
%      basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy2_2013-07\Subject2\Session1\HAAO_STD\Shimmer\';
%      headerPath = [basePath 'SenyssorData_Waist.header'];
%      dataPath = [basePath 'SensorData_Waist.data'];
%      ekfData = ekfDataHandle(headerPath, dataPath);
%      ekfData = ekfData.load; 

properties (SetAccess = public)
%     properties % the following are inheritated from ADataHandle
%         % handle information
%         filepathHeader = '';
%         filepathData = '';
%     end

    %Pointer to the exersise handle this data belongs to so we can trace
    %back up the chain if needed
    exercise_hndl = [];


    % header 
    date = []; % date that the data collection occured
    name = []; % name (or position) of the shimmer
    sensorId = []; % for shimmers, this would be the bluetooth id
    samplingRate = []; % sampling rate of the IMU
    sensorType = []; % the type of sensor used
    accelerometerCalibration = [];
    gyroscopeCalibration = [];

    % data 
    data = []; % data as is from the data file

    % data (reparsed from the data array)
    systemTime = [];
    time = []; 
    dt = [];
    accelerometerUncalibrated = [];
    accelerometerCalibrated = [];
    gyroscopeUncalibrated = [];
    gyroscopeCalibrated = [];

    % data from the file structure that is currently not being stored:
    % FrameSeqNumCalibrated
    % TimeStampUncalibrated,TimeStampCalibrated,SystemMsTimeStampCalibrated,SystemTimeStampCalibrated,
end % properties

methods
    function obj = imuDataHandle(varargin)
        % constructor function. parse the incoming variables
        if nargin > 0
            obj.filepathHeader = varargin{1};
            obj.filepathData = varargin{2};
        end
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

        h1 = subplot(2, 1, 1);
        if ~isempty(obj.accelerometerCalibrated)
            plot(time, obj.accelerometerCalibrated);
            ylabel('Linear acceleration [m/s^2]');
        elseif ~isempty(obj.accelerometerUncalibrated)
            plot(time, obj.accelerometerUncalibrated);
            ylabel('Linear acceleration ADC Values');
        end

        h2 = subplot(2, 1, 2);
        if ~isempty(obj.gyroscopeCalibrated)
            plot(time, obj.gyroscopeCalibrated);
            ylabel('Angular velocity [rad/s]');
        elseif ~isempty(obj.gyroscopeUncalibrated)
            plot(time, obj.gyroscopeUncalibrated);
            ylabel('Linear acceleration ADC Values');
        end
        xlabel('Time [s]');

        linkaxes([h1, h2], 'x');
    end
    
    function rotation = rotate2Grav(obj,grav_axis,start_index,modify,cropStart,cropEnd)
        %Finds a window of points with smallest variance and uses it to
        %calculate a rotation such that T*data = norm(data) in grav_axis. 
        %Data is 3xN, grav axis is 'X' 'Y' 'Z' '-X' '-Y' '-Z'
        %If MODIFY is set to 1 it will rotate the acceleration and gyro
        %data for this sensor

        if nargin == 2
           start_index = 1; 
        end
        
        if nargin == 3
            modify = false;
        end
        
        if exist('cropStart', 'var')
            [~, start_index] = findClosestValue(cropStart, obj.time);
            [~, end_index] = findClosestValue(cropEnd, obj.time);
        else
            end_index = length(obj.time);
        end
        
        dat = obj.accelerometerCalibrated(start_index:end_index,:)';
        
            switch grav_axis
                case 'X'
                    fnct = @rot_x_z_0x;
                case 'Y'
                    fnct = @rot_x_z_0y;
                case 'Z'
                    fnct = @rot_x_z_0z;
                otherwise
                    warning('0 acceleration axis not specified, assuming X');
                    fnct = @rot_x_z_0x;
                    grav_axis = 'X';
            end

            beta = nlinfit(dat',zeros(size(dat,2),1),fnct,[0,0]);

            switch grav_axis
                case 'X'
                    rotation = roty(beta(1))*rotz(beta(2));
                case 'Y'
                    rotation = rotx(beta(1))*rotz(beta(2));
                case 'Z'
                    rotation = rotx(beta(1))*roty(beta(2));
                otherwise
                    rotation = eye(3);
            end            
            
            if exist('modify','var') 
               if modify
                   applyRotation(obj, rotation);
               end
            end
            
            function x = rot_x_z_0x(beta,x)
                x = roty(beta(1))*rotz(beta(2))*x';
                x = x(1,:)';
            end
            
            function y = rot_x_z_0y(beta,x)
                y = rotx(beta(1))*rotz(beta(2))*x';
                y = y(2,:)';
            end
            
            function z = rot_x_z_0z(beta,x)
                z = rotx(beta(1))*roty(beta(2))*x';
                z = z(3,:)';
            end
    end
    
    function applyRotation(obj, rotation)
        dat = rotation*obj.accelerometerCalibrated';
        obj.accelerometerCalibrated = dat';
        
        obj.gyroscopeCalibrated = (rotation*(obj.gyroscopeCalibrated'))';
    end
       
    function offset = gyroOffset(obj,window_size,modify)
            %Calculates the gyro offset, if modify is 1 it will
            %subtract it from the data.

            dat = obj.gyroscopeCalibrated';

            best_var = inf;
            best_index = 1;
            N = 1000;
            if N > (size(dat,2)-window_size)
               N = size(dat,2)-window_size;
            end
            for i=1:N

                variance = mean([var(dat(1,i:i+window_size))...
                    var(dat(2,i:i+window_size))...
                    var(dat(3,i:i+window_size))]);

                if(variance < best_var)
                    best_var = variance;
                    best_index = i;
                end
            end

            offset = mean(dat(:,best_index:best_index+window_size),2);

            if exist('modify','var')
               if modify
                  dat = dat - repmat(offset,1,numel(dat(1,:))); 
                  obj.gyroscopeCalibrated = dat';
               end
            end

        end
    
    function cropData(obj, timeStart, timeEnd)
        [timeStartVal, timeStartInd] = findClosestValue(timeStart, obj.time);
        [timeEndVal, timeEndInd] = findClosestValue(timeEnd, obj.time);
        
        % apply the cropping
        obj.time = obj.time(timeStartInd:timeEndInd);
        obj.accelerometerUncalibrated = obj.accelerometerUncalibrated(timeStartInd:timeEndInd,:);
        obj.accelerometerCalibrated = obj.accelerometerCalibrated(timeStartInd:timeEndInd,:);
        obj.gyroscopeUncalibrated = obj.gyroscopeUncalibrated(timeStartInd:timeEndInd,:);
        obj.gyroscopeCalibrated = obj.gyroscopeCalibrated(timeStartInd:timeEndInd,:);
    end
    
    function interpData(obj,time)
       %Interpolates the calibrated data to the given time using 3d splines 
       for i=1:numel(obj)
           hndl = obj(i);
           [old_time,ia,ic] = unique(hndl.time);
           accel = hndl.accelerometerCalibrated(ia,:)';
           gyro = hndl.gyroscopeCalibrated(ia,:)';
           hndl.accelerometerCalibrated = spline(old_time',accel,time)';
           hndl.gyroscopeCalibrated = spline(old_time',gyro,time)';
           hndl.time = time;
       end 
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
                case 'date'
                    if length(xmlValue) == 10
                        % only has yyyy/mm/dd
                        readString = 'yyyy_mm_dd';
                    else
                        readString = 'yyyy_mm_dd_HH_MM_SS';
                    end
                    
                    obj.date = datenum(xmlValue, readString);
                    
                case 'name'
                    obj.name = xmlValue;
                    
                case 'btId'
                    obj.sensorId = xmlValue;
                    
                case 'samplingRate'
                    obj.samplingRate = str2num(xmlValue);
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
                    
                case 'Calibration'
                    % read calibration data
                    for ind_calibration = 1:length(xmlChild)
                        calibStruct = readCalibrationFromParseXML(xmlChild(ind_calibration).Children);
                        
                        switch xmlChild(ind_calibration).Name
                            case 'ACCEL'
                                obj.accelerometerCalibration = calibStruct;
                            case 'GYRO'
                                obj.gyroscopeCalibration = calibStruct;
                            case 'MAG'
                                obj.magnetometerCalibration = calibStruct;
                        end
                    end
            end
        end
    end
    


    function obj = loadData(obj)
        % load the data via parseCSV
        shimmerData = parseCSV(obj.filepathData);
        
        % check the time vector
        day2minMultiplier = 24*60*60;
        %     assessTimeMS = datenum('1970-1-1') + ekfData.SystemMSTimeStamp(1)/(1000*day2minMultiplier);
        assessTimeS = datenum('1970-1-1') + shimmerData.SystemMsTimeStampCalibrated(1)/(day2minMultiplier);
        
        if assessTimeS > now
            % if this statement is true, then it suggests that the timestamp is
            % in milliseconds and should be scaled down by 1000
            shimmerData.SystemSTimeStamp = shimmerData.SystemMsTimeStampCalibrated/1000;
        else
            shimmerData.SystemSTimeStamp = shimmerData.SystemMsTimeStampCalibrated;
        end
        
        % saving the proper loaded data to the obj
        obj.systemTime = shimmerData.SystemMsTimeStampCalibrated;
        obj.time = shimmerData.SystemSTimeStamp(1) + (shimmerData.TimeStampCalibrated - shimmerData.TimeStampCalibrated(1))/1000;
        
        %Fix time in case it loops over
        %Fix timestamp rollover
        for i =2:numel(obj.time)
           if obj.time(i) < obj.time(i-1)
              obj.time(i:end) = obj.time(i:end)-obj.time(i)+obj.time(i-1); 
           end
        end
        
        
        obj.accelerometerUncalibrated = [shimmerData.AccelerometerXUncalibrated shimmerData.AccelerometerYUncalibrated shimmerData.AccelerometerZUncalibrated];
        
        
        
        if isfield(shimmerData,'AccelerometerXCalibrated')
            obj.accelerometerCalibrated = [shimmerData.AccelerometerXCalibrated shimmerData.AccelerometerYCalibrated shimmerData.AccelerometerZCalibrated];
        end
        
        obj.gyroscopeUncalibrated = [shimmerData.GyroscopeXUncalibrated shimmerData.GyroscopeYUncalibrated shimmerData.GyroscopeZUncalibrated];
        if isfield(shimmerData,'GyroscopeXCalibrated')
            obj.gyroscopeCalibrated = [shimmerData.GyroscopeXCalibrated shimmerData.GyroscopeYCalibrated shimmerData.GyroscopeZCalibrated];
        end
        
        obj.dt = mean(diff(obj.time));
    end

    function writeHeader(obj)

    end

    function writeData(obj)

    end
end % methods (private)
end % classdef