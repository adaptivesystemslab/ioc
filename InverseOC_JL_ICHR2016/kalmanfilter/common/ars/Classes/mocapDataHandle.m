classdef mocapDataHandle < ADataHandle
%  mocapDataHandle is a model to handle data parsing from MoCap data files. 
%  Pass in the header and the data file path to the MoCap data file and it 
%  will load both the header and the data structure 
%
%  Usage:
%    mocapData = mocapDataHandle(filepathHeader, filepathData) 
%      filepathHeader - header filepath
%      filepathData - data filepath
%
%  Example:
%    Loading data from file
%      basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy1_2011-11\Subject5\Session1\KEFO_SIT_SLO1\Cortex\';
%      headerPath = [basePath 'MocapData_Cortex.header'];
%      dataPath = [basePath 'MocapData_Cortex.data'];
%      mocapData = mocapDataHandle(headerPath, dataPath);
%      mocapData = mocapData.load; 
    
properties (SetAccess = public)
%     properties % the following are inheritated from ADataHandle
%         % handle information
%         filepathHeader = '';
%         filepathData = '';
%     end

    exercise_hndl = [];

    offTrackVal = 9999999; % if a marker disappears, it is replaced by this val
    datasetName = [];
    
    % header 
    mocapSource = [];
    headerList = []; % marker list that is in the header for the data
    combineArrayList = []; % names for the knee/ankle markers that should be averaged together
    
    date = []; % date that the EKF was ran
    sensorType = []; % type of mocap system used
    samplingRate = []; % data collection sampling rate
    markerNames = {}; % the name of markers used 

    % data 
    data = []; % data as is from the data file

    % data (reparsed from the data array)
    time = []; 
    dt = [];
end % properties

methods
    function obj = mocapDataHandle(varargin)
        % constructor function. parse the incoming variables
        obj.filepathHeader = varargin{1};
        obj.filepathData = varargin{2};
    end

    function obj = load(obj, varargin)
        % load the data from file based on header and data properties
        if nargin > 1
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
        plot(time, obj.Q);
        ylabel('Joint angle [rad]');

        h2 = subplot(2, 1, 2);
        plot(time, obj.dQ);
        ylabel('Joint velocity [rad/s]');

        xlabel('Time [s]');

        linkaxes([h1, h2], 'x');
    end
    
    function cropData(obj, timeStart, timeEnd)
        [timeStartVal, timeStartInd] = findClosestValue(timeStart, obj.time);
        [timeEndVal, timeEndInd] = findClosestValue(timeEnd, obj.time);
        
        % apply the cropping
        obj.time = obj.time(timeStartInd:timeEndInd);
        
        fields = fieldnames(obj.data);
        for i = 1:numel(fields)
            obj.data.(fields{i}) = obj.data.(fields{i})(timeStartInd:timeEndInd, :);
        end
    end
end % methods

methods (Access = private)
    function obj = loadHeader(obj)

        % todo the following block probably should be integrated into the 
        % rest of the header information, but the list of markers that are
        % currently in the headers are informative and does not actually 
        % match the names in the .data file. Instead, we will guess the
        % datasetName based on the filename
        
        [pathPartsStr, pathPartsName, pathPartsExt] = fileparts(obj.filepathHeader);
        strSplitDirPathExercise = strsplit(obj.filepathHeader, filesep);
        
        obj.datasetName = strSplitDirPathExercise{end-5}; % determine the dataset type based on the file path
   
        switch lower(obj.datasetName)
            case {'healthy1', 'lowerbody_healthy1_2011-11'}
                obj.headerList = {'R_Shoulder', 'L_Shoulder', ...
                    'R_ASIS', 'L_ASIS', ...
                    'R_Knee_Med', 'R_Knee_Lat', 'R_Ankle_Med', 'R_Ankle_Lat', ...
                    'L_Knee_Med', 'L_Knee_Lat', 'L_Ankle_Med', 'L_Ankle_Lat', ...
                    'R_Toe', 'R_Heel', ...
                    'L_Toe', 'L_Heel', ...
                    'Hip_Shimmer', 'R_Knee_Shimmer', 'R_Ankle_Shimmer','L_Knee_Shimmer', 'L_Ankle_Shimmer'};
                obj.combineArrayList = {'R_Knee_Medial', 'R_Knee_Lateral', 'R_Ankle_Medial', 'R_Ankle_Lateral'};
                
            case {'healthy2', 'lowerbody_healthy2_2013-07'}
                obj.headerList = {'R_Shoulder', 'L_Shoulder', ...
                    'R_ASIS', 'L_ASIS', ...
                    'R_Knee_Med', 'R_Knee_Lat', 'R_Ankle_Med', 'R_Ankle_Lat', ...
                    'L_Knee_Med', 'L_Knee_Lat', 'L_Ankle_Med', 'L_Ankle_Lat', ...
                    'R_Toe', 'R_Heel', ...
                    'L_Toe', 'L_Heel', ...
                    'Hip_Shimmer', 'R_Knee_Shimmer', 'R_Ankle_Shimmer','L_Knee_Shimmer', 'L_Ankle_Shimmer'};
                obj.combineArrayList = {'R_Knee_Medial', 'R_Knee_Lateral', 'R_Ankle_Medial', 'R_Ankle_Lateral'};
                
            otherwise
        end
        
        % parse the XML file into a structure
        xDoc = parseXML(obj.filepathHeader);

        % pull the proper data from the base layer
        for ind_attributes = 1:length(xDoc.Attributes)
            xmlName = xDoc.Attributes(ind_attributes).Name;
            xmlValue = xDoc.Attributes(ind_attributes).Value;

            switch xmlName
                case 'name'
                    obj.mocapSource = xmlValue;
                    
                case 'date'
                    obj.date = datenum(xmlValue);
                    
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

                case 'Markers'
                    % discard for now. the markers in the current headers
                    % are informative and does not match the actual .data
                    % list
            end
        end
    end

    function obj = loadData(obj)
        % load the data via parseCSV
        cortexDataTemp = parseCSV(obj.filepathData);

        cortexData.Frame = cortexDataTemp.Frame;
        %cortexData.Frame = [];
        
        % check the time vector
        day2minMultiplier = 24*60*60;
        %     assessTimeMS = datenum('1970-1-1') + ekfData.SystemMSTimeStamp(1)/(1000*day2minMultiplier);
        cortexData.SystemMSTimeStamp = cortexDataTemp.SystemMSTimeStamp;
        assessTimeS = datenum('1970-1-1') + cortexData.SystemMSTimeStamp(1)/(day2minMultiplier);
        
        % check for end-of-day wraparound. ie the clock went from 11.59.59 -> 00.00.00
        diffMS = diff(cortexData.SystemMSTimeStamp);
        [minVal, minInd] = min(diffMS);
        
        if assessTimeS > now
            % if this statement is true, then it suggests that the timestamp is
            % in milliseconds and should be scaled down by 1000
            if minVal < 0 % negative diff...
                cortexData.SystemMSTimeStamp(minInd+1:end) = ...
                    cortexData.SystemMSTimeStamp(minInd+1:end) + 24*60*60*1000;
            end
            
            cortexData.SystemSTimeStamp = cortexData.SystemMSTimeStamp/1000;
        else
            if minVal < 0 % negative diff...
                cortexData.SystemMSTimeStamp(minInd+1:end) = ...
                    cortexData.SystemMSTimeStamp(minInd+1:end) + 24*60*60;
            end
            
            cortexData.SystemSTimeStamp = cortexData.SystemMSTimeStamp;
        end
        
        % look for marker skip and interpolate to fill in
        fieldNames = fieldnames(cortexDataTemp);
        
        %Since healthy1 and healthy2 have inconsistent names here we make
        %them consistent
        
        
        
        lengthData = size(cortexData.SystemSTimeStamp, 1);
        
        % uniform the header information
        prefixString = ''; % pull out any information that prefixes the headers
        %     for i = 1:length(fieldNames)
        %         [matchstart,matchend,tokenindices,matchstring,tokenstring,tokenname,splitstring] ...
        %             = regexp(fieldNames{i}, headerList{1});
        %
        %         if ~isempty(matchstart)
        %             prefixString = splitstring{1};
        %             break
        %         end
        %     end
        
        %We rename some fields here so everything is consistent between
        %healthy 1 and healthy2
        for j = 1:length(fieldNames)
            name = fieldNames{j};
           if(~isempty(strfind(lower(name),'_s_knee')))
              index = strfind(lower(name),'_s_knee');
              cortexDataTemp.([name(1:index) 'Knee_Shimmer' name(index+length('_s_knee'):end)]) = ...
                  cortexDataTemp.(name);
              cortexDataTemp = rmfield(cortexDataTemp,name);
           elseif(~isempty(strfind(lower(name),'_s_waist')))
              index = strfind(lower(name),'_s_waist');
              cortexDataTemp.([name(1:index) 'Hip_Shimmer' name(index+length('_s_waist'):end)]) = ...
                  cortexDataTemp.(name);
              cortexDataTemp = rmfield(cortexDataTemp,name);
           elseif(~isempty(strfind(lower(name),'_s_ankle')))
              index = strfind(lower(name),'_s_ankle');
              cortexDataTemp.([name(1:index) 'Ankle_Shimmer' name(index+length('_s_ankle'):end)]) = ...
                  cortexDataTemp.(name);
              cortexDataTemp = rmfield(cortexDataTemp,name);
           elseif(~isempty(strfind(lower(name),'_lat')) && isempty(strfind(lower(name),'_lateral')))
              index = strfind(lower(name),'_lat');
              cortexDataTemp.([name(1:index) 'Lateral' name(index+4:end)]) = ...
                  cortexDataTemp.(name);
              cortexDataTemp = rmfield(cortexDataTemp,name);
           elseif(~isempty(strfind(lower(name),'_med')) && isempty(strfind(lower(name),'_medial')))
              index = strfind(lower(name),'_med');
              cortexDataTemp.([name(1:index) 'Medial' name(index+4:end)]) = ...
                  cortexDataTemp.(name);
              cortexDataTemp = rmfield(cortexDataTemp,name);
           end
        end
        fieldNames = fieldnames(cortexDataTemp);
        
        for i = 1:length(obj.headerList)
            prefixString = '';
            headerString = '';
            
            for j = 1:length(fieldNames)
                [matchstart,matchend,tokenindices,matchstring,tokenstring,tokenname,splitstring] ...
                    = regexp(fieldNames{j}, obj.headerList{i}, 'ignorecase');
                
                if ~isempty(matchstart)
                    prefixString = splitstring{1};
                    headerString = matchstring{1};
                    
                    if ~isempty(headerString)
                        cortexData.([obj.headerList{i} splitstring{2}]) = cortexDataTemp.([prefixString headerString splitstring{2}]);
                    end
                end
            end
        end
        
        fieldNames = fieldnames(cortexData);
        maxStartIndError = 1; % track the points where the markers disappear, to note where to crop
        minEndIndError = lengthData; 
        %markersToCheck = 'L_ASIS|R_ASIS|Knee|Ankle|Shimmer';
%         markersToCheck = 'L_ASIS|R_ASIS|Knee|Ankle';
        markersToCheck = '[A-z]';
        for i = 4:length(fieldNames) % looking at each dof independently
            % however, we only want to match certain markers, since we
            % don't care if the less important markers (such as heel),
            % which tends to occlude easily, disappears
            [matchstart,matchend,tokenindices,matchstring,tokenstring,tokenname,splitstring] ...
                = regexp(fieldNames{i}, markersToCheck);
            
            if isempty(matchstart)
                continue % skip the markers that doesn't match our list of the things we want to check for
            end
            
            % perform the interpolation to fill in gaps in the data
            %         eval(['sensorDataInsp = cortexData.' fieldNames{i} ';']);
            sensorDataInsp = cortexData.(fieldNames{i});
            [clusterArray, startIndError, endIndError] = obj.interpolateCluster(sensorDataInsp, obj.offTrackVal);
            
            for j = 1:size(clusterArray, 1)
                startInd = clusterArray(j, 1) - 1;
                endInd = clusterArray(j, 2) + 1;
                timeFullSeg = cortexData.SystemSTimeStamp(startInd:endInd); % full time array wanted
                timeSeg = cortexData.SystemSTimeStamp([startInd endInd]); % two time points that correspond to edge points
                dataSeg = sensorDataInsp([startInd endInd]); % corresponding data points
                
                % account for cases where there are only 2 points, since
                % splining doesn't seem to work well with them
                if length(timeSeg) == 2 && length(timeFullSeg) == 3
                    %                 dataSplined = [dataSeg(1) mean(dataSeg) dataSeg(2)];
                    dataSplined = interp1(timeSeg, dataSeg, timeFullSeg, 'linear');
                else
                    dataSplined = spline(timeSeg, dataSeg, timeFullSeg);
                end
                
                cortexData.(fieldNames{i})(startInd:endInd) = dataSplined;
                %             eval(['cortexData.' fieldNames{i} '(startInd:endInd) = dataSplined;']);
            end
            
            if maxStartIndError < startIndError
                maxStartIndError = startIndError;
            end
            
            if minEndIndError > endIndError && endIndError > 1 % endIndError > 1
                minEndIndError = endIndError - 1;
            end
        end
        
% % %         % commit the crop to remove edge data that does not have proper marker
% % %         % data, using maxStartIndError and minEndIndError, then  consolidate
% % %         % the cortex data
% % %         maxStartIndError = maxStartIndError + 1; % remove the first and last points anyway
% % %         minEndIndError = minEndIndError - 1;
% % %         for i = 1:length(fieldNames)
% % %             cortexData.(fieldNames{i}) = cortexData.(fieldNames{i})(maxStartIndError:minEndIndError);
% % %         end
        
        % combine the medial and lateral joints to create an average location between the two markers to get joint centre
        % combine the medial and lateral joints
        %lengthData = size(cortexData.SystemSTimeStamp, 1);
        %cortexData.R_Knee = mean(reshape([cortexData.(obj.combineArrayList{1}) cortexData.(obj.combineArrayList{2})], lengthData, 3, 2), 3);
        %cortexData.R_Ankle = mean(reshape([cortexData.(obj.combineArrayList{3}) cortexData.(obj.combineArrayList{4})], lengthData, 3, 2), 3);
 
        % reparse them explicitly
        obj.data = cortexData;
        obj.time = cortexData.SystemSTimeStamp;
    end

    function writeHeader(obj)

    end

    function writeData(obj)

    end
end % methods (private)
methods(Static)
function [clusterArray, startIndError, endIndError] = interpolateCluster(sensorDataCol, offTrackVal)
    % find the number of "time clusters" where the markers become
    % untrackable, and note their starting and ending points. These points
    % will get interpolated
    
    % clusterArray[1] = start
    % clusterArray[2] = end
    
    clusterInd = 0;
    clusterArray = [];
    
    offTrackArray = find(sensorDataCol == offTrackVal);
    for i = 1:length(offTrackArray)
        if isempty(clusterArray)
            % first entry
            clusterInd = clusterInd + 1;
            
            clusterArray(clusterInd, 1) = offTrackArray(i); 
            clusterArray(clusterInd, 2) = offTrackArray(i);             
        else
            if offTrackArray(i) == clusterArray(clusterInd, 2) + 1
                % if next value is incriment from previous
                clusterArray(clusterInd, 2) = offTrackArray(i);
            else
                % new column
                clusterInd = clusterInd + 1;
                
                clusterArray(clusterInd, 1) = offTrackArray(i);
                clusterArray(clusterInd, 2) = offTrackArray(i);
            end
        end
    end
    
    % the following sets of data should be removed by cropping anyway, but
    % we won't process it at this stage

    % if the last cluster goes all the way to the end, drop that cluster,
    % since it's likely caused by actor going off the capture area, and
    % cannot be interpolated properly anyway
    if ~isempty(clusterArray) && clusterArray(clusterInd, 2) == length(sensorDataCol)
        endIndError = clusterArray(end, 1);
        clusterInd = clusterInd - 1;
        clusterArray = clusterArray(1:clusterInd, 1:2);
    else
        endIndError = 0;
    end
    
    % also, if the first cluster is at the beginning, then we can't do
    % anything about that either, so we'll have to drop that
    if ~isempty(clusterArray) && clusterArray(1, 1) == 1
        startIndError = clusterArray(1, 2);
        clusterArray = clusterArray(2:clusterInd, 1:2);
    else
        startIndError = 0;
    end  
end
end
end % classdef