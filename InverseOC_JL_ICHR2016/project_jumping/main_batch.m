function main_batch(specStruct, run_mode, cost_function_names, outputString, variableFactorsInput, pathToRawData, pathToOutput)
    % preamble loading of background parameters
    addpath(genpath(fullfile('..','Symoro')));
    addpath(genpath(fullfile('..','Model')));
    addpath(genpath(fullfile('..','Common')));
    addpath(genpath(fullfile('..','..','Toolboxes')));
    addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\MarkerSwappingEKF\2018_07_04')); % what folders are correct???
    addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\MarkerSwappingEKF\ekf_2018_05_29')); % what folders are correct???
    
    [cost_function_names, sortInd] = sort(cost_function_names);
    
    for i = 1:length(specStruct.ccost_array)
        specStruct.ccost_array{i} = specStruct.ccost_array{i}(sortInd);
    end
    
%     specStruct.datasetSpecs = datasetSpecs(specStruct.dataset); % had to modify this for jump2D data as well

    switch run_mode
        case 'win'
            specStruct.dataset = lower(specStruct.dataset);
            
        case 'sim'
            specStruct.dataset = 'sim';
    end
    
    runSettings.sharcNet = 0;     % is the function running on the SHARCNET cluster?
    runSettings.plotFig = 1;      % should we plot the output figures and save them?
    runSettings.saveResults = 1;  % do we need to keep a copy of the trained classifier?
    runSettings.verbose = 1;      % should we be verbose to the console or not?
    
%     fileStackTemp = loadPatientFilepaths([pathToRawData '\' specStruct.datasetSpecs.dataPathSuffix], specStruct);
    fileStackTemp = {'tmp'};
    
    nowTimeStr = datestr(now, 'yyyymmdd');
    outputInstancePath = [nowTimeStr '_' run_mode(1) num2str(specStruct.patient(1))];
    
    % Create array of jumps to process, row 1 = partNum, row 2 = targNum (/3), row 3 = jumpNum (/12)
%     partList = strsplit(num2str(specStruct.patient),' ');
%     for i_part = 1:length(partList) % append zeros to make all entries two digits, matching file names
%         if(length(partList{i_part})==1)
%             partList{i_part} = ['0' partList{i_part}];
%         end
%     end
%     
%     
%   
%     jumpDataList = [20; % participant (2-22)
%                     2; % target (1-3)
%                     4]; % jump (1-12)
    jumpDataList = [8, 11, 20, 21; % participant (2-22)
                    2,  2,  2, 2; % target (1-3)
                    4,  4,  4, 8]; % jump (1-12)
%     jumpDataList = [2:22; % participant (2-22)
%                     3*ones(1,21); % target (1-3)
%                     1*ones(1,21)]; % jump (1-12)
    
    for ind_fileStack = 1:size(jumpDataList,2) % number of jumps
        currFilestack.tmp = fileStackTemp; %fileStackTemp{ind_fileStack};
%         currFilePath = currFilestack.filePath;
%         currInstName = ['S' num2str(currFilestack.subjectNumber) ...
%             '_S' num2str(currFilestack.sessionNumber) ...
%             '_' currFilestack.exerciseName];
        
        partNum = num2str(jumpDataList(1,ind_fileStack));
        if(numel(partNum)==1)   partNum = ['0' partNum];    end
        targNum = jumpDataList(2,ind_fileStack);
        jumpNum = jumpDataList(3,ind_fileStack);
        currInstName = ['P' partNum '_T' num2str(targNum) '_J' num2str(jumpNum)];
        
%         fprintf('(%u/%u): %s\n', ind_fileStack, length(fileStackTemp), currInstName);
        disp(['P' partNum ', target ' num2str(targNum) ', jump ' num2str(jumpNum)]);
        
        outputPath = fullfile(pathToOutput, specStruct.dataset, outputString, outputInstancePath, currInstName);
        checkMkdir(outputPath);
        
        currFilestack.dataset = specStruct.dataset;
        
        switch currFilestack.dataset
            case 'sim'
                manSeg = 'Segmentation_manual_JL';
                
                jointAngleFile = searchForFileByExt(currFilePath, 'Kinematics_meas*.mat');
                filesToLoad{1} = fullfile(currFilePath, jointAngleFile);
                
                dynamicsFile = searchForFileByExt(currFilePath, 'Dynamics_meas_*.mat');
                filesToLoad{2} = fullfile(currFilePath, dynamicsFile);
                
            case {'healthy1ankle', 'healthy1hip'}
                manSeg = 'Segmentation_manual';
                filesToLoad{1} = fullfile(currFilePath, 'EKF', '2016_11_14_Planar', 'jointAngles.mat');
                
            case 'squats_urfi_2011'
                manSeg = 'Segmentation_manual_JL';
                
                jointAngleFile = searchForFileByExt(currFilePath, 'Kinematics_meas*.mat');
                filesToLoad{1} = fullfile(currFilePath, jointAngleFile);
                
                dynamicsFile = searchForFileByExt(currFilePath, 'Dynamics_meas_*.mat');
                filesToLoad{2} = fullfile(currFilePath, dynamicsFile);
                

            case {'iit_2017', 'iit_2017_2d'}
                manSeg = '';
                filesToLoad{1} = fullfile(currFilePath, 'JointAngles', 'ekf.mat');
                filesToLoad{2} = fullfile(currFilePath, 'mocap_fp', 'mocap1.trc');
                
            case {'jump2D','jump2d'}
                manSeg = '';
                jumpFilePath = pathToRawData;
                filesToLoad{1} = fullfile(jumpFilePath, ['JA_P', partNum, '_2D.mat']);
                filesToLoad{2} = [num2str(targNum) '_' num2str(jumpNum)];
        end
        
        currFilestack.ccost_array = specStruct.ccost_array;
        currFilestack.cconst_array = specStruct.cconst_array;

        manSegLoadPath = ' '; %fullfile(currFilePath, manSeg, 'SegmentData_Manual_Manual.data');

        runSettings.nowTimeStr = nowTimeStr;
        runSettings.variableFactors = variableFactorsInput;
         
        if 0
            if strcmpi(currInstName(end), '1')
                runSettings.variableFactors.normalizationMethod_cf = 'rang';
                param_test = setup_main(filesToLoad, manSegLoadPath, currFilestack, run_mode, runSettings, [], []);
                norm_coeff_calc = [norm_coeff_calc; param_test.coeff_cf.array];
            else
                fprintf('  Skipping...\n');
            end
            
            if ind_fileStack < length(fileStackTemp) 
                continue
            else
                average_val_use = mean(norm_coeff_calc, 1)
                continue
            end
        end
        
        % load the data to know how much of the file needs parsing
        [param] = setup_main(filesToLoad, manSegLoadPath, currFilestack, 'win', runSettings, [], []); % always set up in 'win' mode not 'sim'

%         continue
        
% % %         % now break it down.
% % %         if variableFactorsInput.dataLoadLength == 0
            startLength = 1;
            endLength = length(param.t_full);
% % %         else
% % %             startLength = min(variableFactorsInput.dataLoadLength);
% % %             endLength = max(variableFactorsInput.dataLoadLength);
% % %             
% % %             if endLength > length(param.t_full)
% % %                 endLength = length(param.t_full);
% % %             end
% % %         end
        
        switch run_mode
            case 'win'
                ind = 0;
                indToUse_window = [];
                skipLength = variableFactorsInput.batchWindowBreakdownFactor - variableFactorsInput.batchWindowBreakdownFactorOverlap;
                for ind_windowCount = startLength:skipLength:endLength
                    ind = ind+1;
                    indToUse_window(ind, 1) = ind_windowCount;
                    indToUse_window(ind, 2) = ind_windowCount + variableFactorsInput.batchWindowBreakdownFactor;
                    
                    if indToUse_window(ind, 2) > endLength
                        %         indToUse_window(ind, 2) = length(t_full);
                        indToUse_window(ind, 2) = endLength; % remove that latest entry
                        break % we're at the end of the line
                    end
                end
                
            case 'sim' % no windows to load
                indToUse_window = [];
                indToUse_window(1, 1) = 1;
                indToUse_window(1, 2) = 1;
        end
        
      
        
% % %         try
            windowCount = size(indToUse_window, 1);
            for ind_windowCount = 1:windowCount
                currInstNameUse = [currInstName '_' num2str(indToUse_window(ind_windowCount, 1)) '_' num2str(indToUse_window(ind_windowCount, 2))];
                variableFactors = variableFactorsInput;
                variableFactors.dataLoadLength = indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2);
                runSettings.variableFactors = variableFactors;
                
                main(filesToLoad, manSegLoadPath, outputPath, currInstNameUse, runSettings, currFilestack, run_mode, cost_function_names);
            end
            
% % %         catch err
% % %             errMessageTmp = regexp(err.message,',','split'); % error messages with commas in it
% % %             errMessageFirstComma = errMessageTmp{1};
% % %             errMessage = [errMessageFirstComma ' - ' err.stack(1).file];
% % %             fprintf('Error START: [%s]\n', errMessage);
% % %             
% % %             for err_ind = 1:length(err.stack)
% % % %                             print the error messages to screen
% % %                 err.stack(err_ind)
% % %             end
% % %             
% % %             fprintf('Error END: [%s]\n', errMessage);
% % %         end
        end
    end
