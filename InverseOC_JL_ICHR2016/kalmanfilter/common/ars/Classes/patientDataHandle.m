classdef patientDataHandle < handle
% exerciseDataHandle aggregiates all the data for a given exercise 
% occurance, and provide file I/O, as well as plotting functions. this
% class functions as a container for the other datasets, and should be
% accessed for all data loading and management
%
% This class can be used as part of a sessionDataHandle, or as an array of
% exerciseDataHandles

properties
    % file I/O
    dirPathRaw = ''; % the path to the base data
    datasetName = ''; % the project that this data is associated with
    patientHeaderPath =''; %Path to the header
    % demographic information
    subjectNumber = 0; % subject number
    gender = '';
    age = [];
    height = [];
    weight = [];
    upperLegLength = [];
    lowerLegLength = [];
    
    % surgical info
    surgicalSide = '';
    surgicalType = '';
    surgicalDate = '';
    surgicalCause = '';
    surgicalComplications = '';
    
    % rehab info
    rehabPT = '';
    rehabAdmissionDate = [];
    rehabDischargeDate = [];
    rehabSessionCount = '';
    priorAssistiveDevices = '';
    
    notes = '';
    
    % available data
    exercises = exerciseDataHandle.empty();
    
    % legacy fields
    subjectName = [];
    
end % properties

methods
    function obj = patientDataHandle(varargin)
        
        if nargin == 1
            %We have only the path, set the header path and move on
            obj.dirPathRaw = varargin{1};
            obj.patientHeaderPath = obj.dirPathRaw;
            obj.loadHeader_patient();
            return;
        end
        constHeader = '.header';
        constData = '.data';
        
        % if only one path passed in, it is assumed that it is the to
        % basepath of the exercise, then the rest is by parameter pairs
        obj.dirPathRaw = varargin{1};
        specStruct = varargin{2};
        options = varargin{3};
        
        if strcmpi(obj.dirPathRaw(end), filesep)
            % there is a slash at the end of the dir pathing. remove it
            obj.dirPathRaw = obj.dirPathRaw(1:end-1);
        end
        
        % based on the file pathing, guess the subject number, session and
        % exercise name
        strSplitDirPathExercise = strsplit(filesep, obj.dirPathRaw);
        
        switch lower(strSplitDirPathExercise{end})
            case 'lowerbody_healthy1_2011-11'
                obj.datasetName = 'healthy1';
                
            case 'lowerbody_healthy2_2013-07'
                obj.datasetName = 'healthy2';
                
            case 'lowerbody_tri1_2012-10'
                obj.datasetName = 'tri';
        end        
        
        % parse the patient header file
        patientHeaderPath = fullfile(obj.dirPathRaw, ['Subject' num2str(specStruct.patient)], ...
            ['Subject' num2str(specStruct.patient) '.header']);
        obj.patientHeaderPath = patientHeaderPath;
        obj.loadHeader_patient();
        
        % populate exercise based on specStruct
        fileStack = loadPatientFilepaths(obj.dirPathRaw, specStruct);
        
        % load data
        obj.exercises = exerciseDataHandle.empty(numel(fileStack),0);
        for ind_fileStack = 1:length(fileStack)
            %Load the exersise data
            exercise_hndl = exerciseDataHandle(fileStack{ind_fileStack}.filePath, options);
            %Point to this patient in exersise so we can trace back if
            %needed to patient information
            exercise_hndl.patient_hndl = obj;
            obj.exercises(ind_fileStack) = exercise_hndl;
        end
        
        % legacy fields
        obj.subjectName = obj.subjectNumber;
    end
    
    function obj = load(obj,varargin)
        %Works even if we have an array of patients
        for i=1:numel(obj)
            for ind_exercises = 1:length(obj(i).exercises)
                %Load the exersise data
                obj(i).exercises(ind_exercises).load;
            end
        end
    end
    
    function obj = loadHeader_patient(obj)
        % parse the XML file into a structure
        xDoc = parseXML(obj.patientHeaderPath);
        
        % pull the proper data from the base layer
        for ind_attributes = 1:length(xDoc.Attributes)
            xmlName = xDoc.Attributes(ind_attributes).Name;
            xmlValue = xDoc.Attributes(ind_attributes).Value;
            
            switch xmlName
                case 'age'
                    obj.age = str2num(xmlValue);
                    
                case 'gender'
                    obj.gender = upper(xmlValue);
                    
                case 'id'
                    
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
                    
                case 'Anthropometrics'
                    for ind_childAttributes = 1:length(xmlAttributes)
                        xmlChildName = xmlAttributes(ind_childAttributes).Name;
                        xmlChildValue = xmlAttributes(ind_childAttributes).Value;
                        
                        switch xmlChildName
                            case 'height'
                                obj.height = str2num(xmlChildValue);
                            case 'weight'
                                obj.weight = str2num(xmlChildValue);
                            case 'upperLegLength'
                                obj.upperLegLength = str2num(xmlChildValue);
                            case 'lowerLegLength'
                                obj.lowerLegLength = str2num(xmlChildValue);
                        end
                    end
                    
                case 'Surgery'
                    for ind_childAttributes = 1:length(xmlAttributes)
                        xmlChildName = xmlAttributes(ind_childAttributes).Name;
                        xmlChildValue = xmlAttributes(ind_childAttributes).Value;
                        
                        switch xmlChildName
                            case 'type'
                                obj.surgicalType = xmlChildValue;
                            case 'side'
                                obj.surgicalSide = xmlChildValue;
                            case 'date'
                                obj.surgicalDate = datenum(xmlChildValue, 'yyyy_mm_dd');
                            case 'cause'
                                obj.surgicalCause = xmlChildValue;
                            case 'complications'
                                obj.surgicalComplications = xmlChildValue;
                        end
                    end
                    
                case 'Rehabilitation'
                    for ind_childAttributes = 1:length(xmlAttributes)
                        xmlChildName = xmlAttributes(ind_childAttributes).Name;
                        xmlChildValue = xmlAttributes(ind_childAttributes).Value;
                        
                        switch xmlChildName
                            case 'PT'
                                obj.rehabPT = str2num(xmlChildValue);
                            case 'admissionDate'
                                obj.rehabAdmissionDate = datenum(xmlChildValue, 'yyyy_mm_dd');
                            case 'priorAssistiveDevices'
                                obj.priorAssistiveDevices = xmlChildValue;
                            case 'dischargeDate'
                                if isempty(xmlChildValue)
                                   xmlChildValue = '2014_01_01'; 
                                end
                                obj.rehabDischargeDate = datenum(xmlChildValue, 'yyyy_mm_dd');
                            case 'session'
                                obj.rehabSessionCount = str2num(xmlChildValue);
                        end
                    end
                        
                case 'Notes'
                    for ind_childAttributes = 1:length(xmlAttributes)
                        xmlChildName = xmlAttributes(ind_childAttributes).Name;
                        xmlChildValue = xmlAttributes(ind_childAttributes).Value;
                        
                        switch xmlChildName
                            case 'text'
                                obj.notes = xmlChildValue;
                        end
                    end
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