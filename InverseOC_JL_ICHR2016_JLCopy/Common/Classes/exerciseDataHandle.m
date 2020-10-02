classdef exerciseDataHandle < matlab.mixin.Copyable
% exerciseDataHandle aggregiates all the data for a given exercise 
% occurance, and provide file I/O, as well as plotting functions. this
% class functions as a container for the other datasets, and should be
% accessed for all data loading and management
%
% This class can be used as part of a sessionDataHandle, or as an array of
% exerciseDataHandles

properties
    % file I/O
    dirPathExercise = ''; % the path to the base data
    datasetName = ''; % the project that this data is associated with
    
    %Pointer to the PATIENT object
    patient_hndl = 0;
    
    % demographic information
    subjectNumber = 0; % subject number
    sessionNumber = 0; % session
    exerciseName = ''; % the name of the folder
    
    % exercise info
    exerciseType = ''; % the type of exercise
    exerciseNameLong = '';
    kinematicDirection = ''; % the child data is parsed as a string, so this converts it back to an array
    initPosture = '';
    description = '';
    
    % identifier strings
    subjSessExer = '';
    
    % header files to load
    filepathHeaderMocap = [];
    filepathHeaderImu = [];
    filepathHeaderEkf = [];
    filepathHeaderSegManual = [];
    filepathHeaderSegCrop = [];
    filepathHeaderSegAlg = [];
    
    % available data
    dataMocap = [];
    dataImu = [];
    dataEkf = [];
    dataSegCrop = [];
    dataSegManual = [];
    dataSegAlg = [];
    
    % legacy fields
    subjectName = [];
    session = [];
end % properties

methods
    function obj = exerciseDataHandle(varargin)
        constHeader = '.header';
        constData = '.data';

        if isstruct(varargin{1})
            % if it is a struct, there are existing data we can load into
            % this class. 
            specStruct = varargin{1};
            obj.dirPathExercise = specStruct.filePath;
            
            % if the user passed in a specific path seq, then we probably 
            % don't want to load all the default strings
            populateDefaultStrings = 0; 
        else
            % if only one path passed in, it is assumed that it is the to
            % basepath of the exercise, then the rest is by parameter pairs
            obj.dirPathExercise =  varargin{1};
            specStruct = varargin{2};
            
            % on the other hand, passing in just the base path probably
            % mean they've passed in an options file, which then we should
            % load everything
            populateDefaultStrings = 1;
        end
        
        if strcmpi(obj.dirPathExercise(end), filesep)
            % there is a slash at the end of the dir pathing. remove it
            obj.dirPathExercise = obj.dirPathExercise(1:end-1);
        end
        
        % based on the file pathing, guess the subject number, session and
        % exercise name
        strSplitDirPathExercise = strsplit(obj.dirPathExercise, filesep);
        
        switch lower(strSplitDirPathExercise{end-3})
            case 'lowerbody_healthy1_2011-11'
                subjectDatasetName = 'healthy1';
                
            case 'lowerbody_healthy2_2013-07'
                subjectDatasetName = 'healthy2';
                
            case 'lowerbody_tri1_2012-10'
                subjectDatasetName = 'tri';
                
            case 'lowerbody_stjoseph1_2013-02'
                subjectDatasetName = 'stjoseph1';
                
            case 'grh'
                subjectDatasetName = 'grh';
                
            otherwise
                subjectDatasetName = 'undefined';
        end
        
        subjectId = str2num(strSplitDirPathExercise{end-2}(8:end)); % assumes the word 'subject' precedes it
        subjectSession = str2num(strSplitDirPathExercise{end-1}(8:end)); % assume the word 'session' precedes it
        subjectExerciseName = strSplitDirPathExercise{end};
        subjectExerciseType = subjectExerciseName;
        for ind_exerciseNameLength = 1:length(subjectExerciseName)
            % want to remove the trailing number on the exercise name
            intCheck = str2num(subjectExerciseName(ind_exerciseNameLength:end));
            if ~isempty(intCheck)
                subjectExerciseType = subjectExerciseName(1:ind_exerciseNameLength-1);
                break
            end
        end
        
        % if there is an exerciseDescription folder, load the proper
        % exercise description for this movement
        exerciseDescPath = fullfile(obj.dirPathExercise, '..', '..', '..', 'ExerciseDescription');
        exerciseDescDir = dir(fullfile(exerciseDescPath, [subjectExerciseType(1:8) '.header']));
        if ~isempty(exerciseDescDir) && ~strcmpi(computer, 'GLNXA64')
            % load the contents, but not if we're running on SHARCNET
            loadHeader_exerciseDesc(obj, fullfile(exerciseDescPath, exerciseDescDir.name));
        end
        
        % search in the directory for the default pathing to the individual
        % files. Locate the header files if possible, given default
        % assumptions
        dirPathMocap = fullfile(obj.dirPathExercise, 'Cortex');
        if ~exist(dirPathMocap, 'dir') || ~populateDefaultStrings
            filepathHeaderMocap = [];
        else
            filepathHeaderMocap = obj.extFileSearch(dirPathMocap, constHeader);
        end
        
        dirPathImu = fullfile(obj.dirPathExercise, 'IMU'); % potentially have different path names
        dirPathImu2 = fullfile(obj.dirPathExercise, 'Shimmer');
        if (~exist(dirPathImu, 'dir') && ~exist(dirPathImu2, 'dir')) || ~populateDefaultStrings
            filepathHeaderImu = [];
        else
            % iterate through the shimmer folder to pull out all the header
            % files for all IMU devices
            if exist(dirPathImu, 'dir')
                filepathHeaderImu = obj.extFileSearch(dirPathImu, constHeader);
            elseif exist(dirPathImu2, 'dir')
                filepathHeaderImu = obj.extFileSearch(dirPathImu2, constHeader);
            end
        end
        
        dirPathEkf = fullfile(obj.dirPathExercise, 'EKF');
        if ~exist(dirPathEkf, 'dir') || ~populateDefaultStrings
            filepathHeaderEkf = [];
        else
            % load the latest EKF folder
            ekfSubfolderDir = dir(dirPathEkf);
            dirPathEkf = fullfile(dirPathEkf, ekfSubfolderDir(end).name); % assume the last one in the dir is the one we want
%              dirPathEkf = fullfile(dirPathEkf, ekfSubfolderDir(3).name);
            
            filepathHeaderEkf = obj.extFileSearch(dirPathEkf, constHeader); % this one we'll use data for now
            
            if isempty(filepathHeaderEkf)
                % DEBUG some EKF data doesn't have header data right now,
                % such as healthy1. This checks specifically for that case
                filepathHeaderEkf = obj.extFileSearch(dirPathEkf, constData, 'Joint'); % this one we'll use data for now
            end
        end
        
        dirPathSegManual = fullfile(obj.dirPathExercise, 'Segmentation_manual');
        if ~exist(dirPathSegManual, 'dir') || ~populateDefaultStrings
            filepathHeaderSegManual = [];
        else
            filepathHeaderSegManual = obj.extFileSearch(dirPathSegManual, constHeader);
        end
        
        dirPathSegCrop = fullfile(obj.dirPathExercise, 'Segmentation_cropping');
        if ~exist(dirPathSegCrop, 'dir') || ~populateDefaultStrings
            filepathHeaderSegCrop = [];
        else
            filepathHeaderSegCrop = obj.extFileSearch(dirPathSegCrop, constHeader);
        end
        
        dirPathSegAlg = fullfile(obj.dirPathExercise, 'Segmentation_algorithmic');
        if ~exist(dirPathSegAlg, 'dir') || ~populateDefaultStrings
            filepathHeaderSegAlg = [];
        else
            % load the latest alg seg folder
            algSegSubfolderDir = dir(dirPathSegAlg);
            dirPathSegAlg = fullfile(dirPathSegAlg, algSegSubfolderDir(end).name); % assume the last one in the dir is the one we want
            
            filepathHeaderSegAlg = obj.extFileSearch(dirPathSegAlg, constHeader);
        end
        
        % replace the default values with the actual values if they're
        % being passed in explicitly
        p = inputParser;
        p.KeepUnmatched = true;
 
        addOptional(p, 'datasetName', subjectDatasetName);
        addOptional(p, 'subjectName', subjectId);
        addOptional(p, 'session', subjectSession);
        addOptional(p, 'exerciseName', subjectExerciseName);
        addOptional(p, 'exerciseType', subjectExerciseType);
        
        addOptional(p, 'headerMocap', filepathHeaderMocap);
        addOptional(p, 'headerImu', filepathHeaderImu);
        addOptional(p, 'headerEkf', filepathHeaderEkf);
        addOptional(p, 'headerSegManual', filepathHeaderSegManual); 
        addOptional(p, 'headerSegCrop', filepathHeaderSegCrop);
        addOptional(p, 'headerSegAlg', filepathHeaderSegAlg);

        parse(p, specStruct); % perform the file checking
        
        % save the data passed in, or use default values
        obj.datasetName = p.Results.datasetName;
        obj.subjectNumber = p.Results.subjectName;
        obj.sessionNumber = p.Results.session;
        obj.exerciseName = p.Results.exerciseName;
        obj.exerciseType = p.Results.exerciseType;
        
        obj.filepathHeaderMocap = p.Results.headerMocap;
        obj.filepathHeaderImu = sort(p.Results.headerImu);
        obj.filepathHeaderEkf = p.Results.headerEkf;
        obj.filepathHeaderSegManual = p.Results.headerSegManual;
        obj.filepathHeaderSegCrop = p.Results.headerSegCrop;
        obj.filepathHeaderSegAlg = p.Results.headerSegAlg;        
        
        % ---
        
        % legacy fields
        obj.subjectName = obj.subjectNumber;
        obj.session = obj.sessionNumber;
        
        % populate pre-canned identifier strings
        obj.subjSessExer = ['Subj' num2str(obj.subjectName) '_' ...
            'Sess' num2str(obj.session) '_' ...
            obj.exerciseName];
    end
    
    function obj = load(obj,varargin)
        if ~isempty(obj.filepathHeaderMocap)
            obj.dataMocap = arsLoader(obj.filepathHeaderMocap);
            obj.dataMocap.exercise_hndl = obj;
        end
        
        if ~isempty(obj.filepathHeaderImu)
            % if there is more than 1 imu. the load list is sorted by file
            % name and will load in that sequence
            if length(obj.filepathHeaderImu) > 1
                obj.dataImu = imuDataHandle.empty(numel(obj.filepathHeaderImu),0);
                imuCounter = 0;
                for ind_imu = 1:length(obj.filepathHeaderImu)
                    if exist(obj.filepathHeaderImu{ind_imu}, 'file')
                        imuCounter = imuCounter + 1;
                        obj.dataImu(imuCounter) = arsLoader(obj.filepathHeaderImu{ind_imu});
                        obj.dataImu(imuCounter).exercise_hndl = obj;
                    end
                end
            else
                % just one imu
                obj.dataImu = arsLoader(obj.filepathHeaderImu);
                obj.dataImu.exercise_hndl = obj;
            end
        end
        
        if ~isempty(obj.filepathHeaderEkf)
            obj.dataEkf = arsLoader(obj.filepathHeaderEkf);
        end
        
        if ~isempty(obj.filepathHeaderSegManual)
            obj.dataSegManual = arsLoader(obj.filepathHeaderSegManual);
        end
        
        if ~isempty(obj.filepathHeaderSegCrop)
            obj.dataSegCrop = arsLoader(obj.filepathHeaderSegCrop);
        end
        
        if ~isempty(obj.filepathHeaderSegAlg)
            obj.dataSegAlg = arsLoader(obj.filepathHeaderSegAlg);
        end
        
        %In some cases the patient file specifies operation on left leg but
        %imu names were always Waist RKnee RAnkle so here we rename those
        %caes to LKnee and LAnkle
        
        if isfield(obj.patient_hndl,'surgicalSide') && numel(obj.dataImu) == 3 && strcmp(obj.patient_hndl.surgicalSide,'L')
           
           %Gotta rename the IMUs
           for i=1:numel(obj.dataImu)
               index = strfind(obj.dataImu(i).name,'R');
              if(~isempty(index))
                  obj.dataImu(i).name(index) = 'L';
              end
           end
            
        end
        
        
    end
    
    function h = plot(obj)
        colourCrop = 'k';
        colourSegmentGood = 'g';
        colourSegmentBad = 'r';
        
        h = figure;
        
        timeToPlot = obj.dataEkf.time;
        dataToPlot = obj.dataEkf.Q;
        veloToPlot = obj.dataEkf.dQ;
        
        % pull out the cropping results
        if ~isempty(obj.dataSegCrop)
            cropValStart = obj.dataSegCrop.timeStart;
            cropValEnd = obj.dataSegCrop.timeEnd;
        else
            cropValStart = [];
            cropValEnd = [];
        end
            
        % and segmentation results
        if ~isempty(obj.dataSegManual)
            algValStart = obj.dataSegManual.timeStart;
            algValEnd = obj.dataSegManual.timeEnd;
            useArray = find(obj.dataSegManual.use);
            useNotArray = setxor(1:length(algValStart), useArray);
        else
            algValStart = [];
            algValEnd = [];
        end
        
        % plotting joint angles
        h1 = subplot(2, 1, 1);
        plot(timeToPlot, dataToPlot);
        
        hold on
        
        h_t = title([obj.datasetName '_subj' num2str(obj.subjectName) ...
                '_sess' num2str(obj.session) ...
                '_' obj.exerciseName]);
        set(h_t, 'Interpreter', 'none')

        ylabel('Joint angles [rad]');    
        
        % and cropping/segmentation
        if ~isempty(cropValStart)
            plotBoxes(h, cropValStart, cropValEnd, colourCrop, 0);
        end
        
        if ~isempty(algValStart)
            plotBoxes(h, algValStart(useArray), algValEnd(useArray), colourSegmentGood, 0.01);
            plotBoxes(h, algValStart(useNotArray), algValEnd(useNotArray), colourSegmentBad, 0.01);
        end
        
        % now plotting joint velocity
        h2 = subplot(2, 1, 2);
        plot(timeToPlot, veloToPlot);
        
        ylabel('Joint velo [rad]');
        
        if ~isempty(cropValStart)
            plotBoxes(h, cropValStart, cropValEnd, colourCrop, 0);
        end
        
        if ~isempty(algValStart)
            plotBoxes(h, algValStart(useArray), algValEnd(useArray), colourSegmentGood, 0.01);
            plotBoxes(h, algValStart(useNotArray), algValEnd(useNotArray), colourSegmentBad, 0.01);
        end
        
        linkaxes([h1, h2], 'x');
        
        xlabel('Time [s]');
    end
    
    function [timeStart, timeEnd] = cropTime(obj,cropType,i)
        if ~exist('cropType', 'var')
            cropType = 'any';
        end
        
        if ~exist('i', 'var')
            i = 1;
        end
        
        if ~isempty(obj(i).dataSegCrop) && (strcmpi(cropType, 'any') || strcmpi(cropType, 'segCrop'))
            timeStart = obj(i).dataSegCrop.timeStart;
            timeEnd = obj(i).dataSegCrop.timeEnd;
            
        elseif ~isempty(obj(i).dataSegManual)
            tol = 5; % crop with a 5 second gap
            timeStart = obj(i).dataSegManual.timeStart(1) - tol;
            timeEnd = obj(i).dataSegManual.timeEnd(end) + tol;
            
        else
            fprintf('No cropping data loaded \n');
            timeStart = [];
            timeEnd = [];
        end
    end
    
    function obj = cropData(obj,cropType)
        % apply the cropping segmentation to each data sequence available
        if ~exist('cropType', 'var')
            cropType = 'any';
        end
        
        for i=1:numel(obj)
            [timeStart, timeEnd] = cropTime(obj,cropType,i);
            if isempty(timeStart)
                fprintf('No cropping data loaded \n');
                return
            end

            if ~isempty(obj(i).dataMocap)
                obj(i).dataMocap.cropData(timeStart, timeEnd);
            end

            if ~isempty(obj(i).dataImu)
                for ind = 1:length(obj(i).dataImu)
                    obj(i).dataImu(ind).cropData(timeStart, timeEnd);
                end
            end

            if ~isempty(obj(i).dataEkf)
                obj(i).dataEkf.cropData(timeStart, timeEnd);
            end

            % no cropping should be applied to the segment data
        end
    end
    
    function obj = loadHeader_exerciseDesc(obj, exerciseDescHeaderPath)
        
        % parse the XML file into a structure
        xDoc = parseXML(exerciseDescHeaderPath);
        
        % pull the proper data from the base layer
        for ind_attributes = 1:length(xDoc.Attributes)
            xmlName = xDoc.Attributes(ind_attributes).Name;
            xmlValue = xDoc.Attributes(ind_attributes).Value;
            
            switch xmlName
                case 'short'
                    
                case 'long'
                    obj.exerciseNameLong = xmlValue;
            end
        end
        
        for ind_attributes = 1:length(xDoc.Children)
            xmlName = xDoc.Children(ind_attributes).Name;
            xmlData = xDoc.Children(ind_attributes).Data;
            xmlChild = xDoc.Children(ind_attributes).Children;
            
            % sometimes the xmlChild will be empty if the attribute itself
            % is empty, causing an unexpected struct situation
            if isempty(xmlChild)
                xmlChild.Data = [];
            end
            
            switch xmlName
                case '#text'
                    if sum(isspace(xmlData)) > 0
                        % parsed a blank. do nothing. not sure why this
                        % happens
                    end
                    
                case 'kinematicDirection'
                    obj.kinematicDirection = xmlChild.Data; % the child data is parsed as a string, so this converts it back to an array
                    
                case 'initPosture'
                    obj.initPosture = xmlChild.Data;
                    
                case 'description'
                    obj.description = xmlChild.Data;
            end
        end
    end
end % methods

methods(Static)
    function fileExtList = extFileSearch(dirPath, extToLocate, prefixToLocate)
        % given a directory path, search its contents for all .[ext] (some
        % extension). The likely usage of this is to search for .header and
        % .data files for file pathing purposes. fileExtList will be
        % returned as a string, unless there is multiple entries of a given
        % file found, in which case it would be returned as a 
        
        if ~exist('prefixToLocate', 'var')
            prefixToLocate = [];
        end
        
        fileExtList = {};
        fileExtCounter = 0;
        dirPathDir = dir(dirPath);
        for ind_dirPath = 1:length(dirPathDir)
            currFile = fullfile(dirPath, dirPathDir(ind_dirPath).name);
            [pathstr, name, ext] = fileparts(currFile);
            
            if strcmpi(extToLocate, ext) && ...
                    (isempty(prefixToLocate) || strcmpi(name(1:length(prefixToLocate)), prefixToLocate))
                fileExtCounter = fileExtCounter + 1;
                fileExtList{fileExtCounter} = currFile;
            end
        end
        
        if fileExtCounter == 1
            fileExtList = fileExtList{1};
        end
    end
end
end % classdef