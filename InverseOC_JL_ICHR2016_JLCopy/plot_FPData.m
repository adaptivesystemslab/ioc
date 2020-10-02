function plot_FPData
    % for each subject, load the FP data and examine the balance over time

%     filename_forces = 'C:\Documents\aslab\data\Squats_TUAT_2011-06\Subject1\Session1\SQUA_STD_NON2\Forces\subject1_02.forces';
%     filename_seg = 'C:\Documents\aslab\data\Squats_TUAT_2011-06\Subject1\Session1\SQUA_STD_NON2\Segmentation_manual_MK\SegmentData_Manual_Manual.data';
%     loadCOPDataAndPlot(filename_forces, filename_seg)

    main_batch;
end

function main_batch(specStruct, outputString)
    % preamble loading of background parameters
    addpath(genpath(fullfile('..','Symoro')));
    addpath(genpath(fullfile('..','Model')));
    addpath(genpath(fullfile('..','Common')));
    addpath(genpath(fullfile('..','..','Toolboxes')));
    
    if ~exist('specStruct', 'var')
%         specStruct.dataset = 'Squats_TUAT_2015';
%         specStruct.patient = [1 2 3 4 5 6 7 8];
%         specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
%         specStruct.exerciseAcceptPrefixInstance = 1;
        
        specStruct.dataset = 'Squats_TUAT_2011';
        specStruct.patient = [1 2 3 4 5 6];
        specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
        specStruct.exerciseAcceptPrefixInstance = 0;
    end
    
    if ~exist('outputString', 'var')
        outputString = 'testing';
    end
    
    run_mode = 'win';
    specStruct.datasetSpecs = datasetSpecs(specStruct.dataset);
    specStruct.dataset = lower(specStruct.dataset);
    
    runSettings.sharcNet = 0;     % is the function running on the SHARCNET cluster?
    runSettings.plotFig = 1;      % should we plot the output figures and save them?
    runSettings.saveResults = 1;  % do we need to keep a copy of the trained classifier?
    runSettings.verbose = 1;      % should we be verbose to the console or not?
    
    pathToRawData = fullfile('C:', 'Documents', 'aslab', 'data', specStruct.datasetSpecs.dataPathSuffix);
    outputBasePath = fullfile('C:', 'Documents', 'MATLABResults', 'IOCProject', specStruct.dataset);
    
    fileStackTemp = loadPatientFilepaths(pathToRawData, specStruct);
    
    nowTimeStr = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
    outputInstancePath = [nowTimeStr '_' run_mode '_' num2str(specStruct.patient(1))];

    for ind_fileStack = 1:length(fileStackTemp)
        currFilestack = fileStackTemp{ind_fileStack};
        currFilePath = currFilestack.filePath;
        currInstName = ['Subj' num2str(currFilestack.subjectNumber) ...
            '_Sess' num2str(currFilestack.sessionNumber) ...
            '_Exer' currFilestack.exerciseName];        

        outputPath = fullfile(outputBasePath, outputString, outputInstancePath, currInstName);
        checkMkdir(outputPath);
        
        currFilestack.dataset = specStruct.dataset;

        switch currFilestack.dataset
            case 'squats_tuat_2011'
                manSeg = 'Segmentation_manual_MK';
                
                jointAngleFile = searchForFileByExt(fullfile(currFilePath, 'JointAngles', 'IK_2016-01_MK'), 'sub*.mat');
                filesToLoad{1} = fullfile(fullfile(currFilePath, 'JointAngles', 'IK_2016-01_MK'), jointAngleFile);
                
                mocapFile = searchForFileByExt(fullfile(currFilePath, '..', 'prm'), 'sub*.mat');
                filesToLoad{2} = fullfile(fullfile(currFilePath, '..', 'prm'), mocapFile);
                
                jointAngleFile = searchForFileByExt(fullfile(currFilePath, 'Forces'), 'sub*.forces');
                forceFile = fullfile(fullfile(currFilePath, 'Forces'), jointAngleFile);
                
                if isempty(jointAngleFile)
                    continue
                end
                
                switch  runSettings.sharcNet
                    case 0 
                        filesToLoad{3} = 'C:\Documents\aslab\data\Squats_TUAT_2015-12\Subject01\Session1\SQUA_STD_NON1\Dynamics_meas_Sb1_Tr1.mat';
            
                    case 1
                        filesToLoad{3} = fullfile('..', '..', 'data', 'Squats_TUAT_2015-12', 'Subject01', 'Session1', 'SQUA_STD_NON1', 'Dynamics_meas_Sb1_Tr1.mat');
                end
                
            case 'squats_tuat_2015'
                manSeg = 'Segmentation_manual_JL';
                
                jointAngleFile = searchForFileByExt(currFilePath, 'Kinematics_meas*.mat');
                filesToLoad{1} = fullfile(currFilePath, jointAngleFile);
                
                dynamicsFile = searchForFileByExt(currFilePath, 'Dynamics_meas_*.mat');
                filesToLoad{2} = fullfile(currFilePath, dynamicsFile);
        end
        
        currFilestack.ccost_array = [];
        currFilestack.cconst_array = [];

        manSegLoadPath = fullfile(currFilePath, manSeg, 'SegmentData_Manual_Manual.data');

        runSettings.nowTimeStr = nowTimeStr;
        runSettings.variableFactors = [];
         
        % load the data to know how much of the file needs parsing
        [param, traj_load, segmentInfo] = setup_main(filesToLoad, manSegLoadPath, currFilestack, run_mode, runSettings, [], []);
        
        fp_data = loadCOPDataAndPlot(forceFile);
        
        h = figure;
        plot(fp_data.t, fp_data.cop_x, 'b', 'DisplayName', 'copx'); hold on
        plot(fp_data.t, fp_data.cop_y, 'r', 'DisplayName', 'copy');
        legend show
        
            % by reps
            if ~isfield(param, 'fatigue_level') || sum(param.fatigue_level)*5 > length(segmentInfo.timeStart)
                % no fatigue level defined, or
                % saved fatigue level is higher than the amount of segment defined
                param.fatigue_level = [0 0 0 0 length(segmentInfo.timeStart)/5];
            end
            entriesToAdd = 1:((param.fatigue_level(1)-1) * 5) + 5;
            if isempty(entriesToAdd)
                entriesToAdd = 0;
            end
            plotGroups = {entriesToAdd};
            for ind_rep = 2:length(param.fatigue_level)
                entriesToAdd = plotGroups{ind_rep-1}(end) + [1:((param.fatigue_level(ind_rep)-1) * 5) + 5];
                if isempty(entriesToAdd)
                    entriesToAdd = 0;
                end
                plotGroups{ind_rep} = entriesToAdd;
            end
        
        colours = {'b', 'g', 'm', 'k', 'c'};
        for ind = 1:length(plotGroups)
            if plotGroups{ind} ~= 0
                plotBoxes(gcf, segmentInfo.timeStart(plotGroups{ind}), segmentInfo.timeEnd(plotGroups{ind}), colours{ind}, 0, 500, 0);
                
                % calc and plot the mean and std for each section
            end
        end
        
        saveas(h, [outputPath '_fp.fig']);
        close(h);
    end
end


function fp_data = loadCOPDataAndPlot(filename_forces)
    dt = 0.01;

    % load force data
    force_struct = readForce(filename_forces);
    
    % calculate COP
%     cop_x = force_struct.data.FX1/force_struct.data.FZ1;
    t = force_struct.data.Sample * dt;
    cop_x = force_struct.data.X1;
    cop_y = force_struct.data.Y1;
    
    fp_data.t = t;
    fp_data.cop_x = cop_x;
    fp_data.cop_y = cop_y;
    
%     parseCSVsettings.skipLastEntry = 0;
%     segData = parseCSV(filename_seg, parseCSVsettings);
%     
%     % load segment data
%     segmentInfo.use = segData.Use;
%     segmentInfo.segmentCount = 1:length(segData.Use);
%     segmentInfo.timeStart = segData.TimeStart * dt;  % convert frames to seconds
%     segmentInfo.timeEnd = segData.TimeEnd * dt;
%     
%     figure;
%     plot(t, cop_x, 'b', 'DisplayName', 'copx'); hold on
%     plot(t, cop_y, 'r', 'DisplayName', 'copy'); 
%     legend show
%     plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k');
end