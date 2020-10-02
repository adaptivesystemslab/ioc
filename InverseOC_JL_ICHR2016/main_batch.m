function main_batch(specStruct, run_mode, cost_function_names, outputString, variableFactorsInput, pathToRawData, pathToOutput)
    % preamble loading of background parameters
       
%     addpath(genpath(fullfile('Common')));
%     addpath(genpath('../../../PoseEstimation/kalmanfilter/ik_framework/common'));
%     addpath(genpath('../../../PoseEstimation/kalmanfilter/ik_framework/instance_expressiveioc'));
%     addpath(genpath('../../../PoseEstimation/kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));
    
%     addpath(genpath(fullfile('Common')));
%     addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\common'));
%     addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc'));
%     addpath(genpath('D:\aslab_gitlab\kalmanfilter\General_FKEKF\DynamicsModelMatlab\MatlabWrapper'));
%     sourceDataTrcFolder = 'D:/aslab/data/expressiveiocData/dataTrc/PamelasData/';
%     sourceDataMatFolder = 'D:/aslab/data/expressiveiocData/dataMat/2019_04_11_rightarm3/matEkfIk';

    addpath(genpath(fullfile('Common')));
    addpath(genpath('../../../kalmanfilter/ik_framework/common'));
    addpath(genpath('../../../kalmanfilter/ik_framework/instance_expressiveioc'));
    addpath(genpath('../../../kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));
    
    sourceDataTrcFolder = '../../expressiveiocData/dataTrc/PamelasData/';
    sourceDataMatFolder = '../../expressiveiocData/dataMat/2019_04_11_rightarm3/matEkfIk';

    fileStackTempTemp = dir(sourceDataMatFolder);
    fileStackTempInd = 0;
    
    calculateNormalization = 0;
    
    segmentsOfInterest = createTaskRepetitionTimeStruct;
    
    norm_coeff_calc = [];
    
    for i = 1:length(fileStackTempTemp)
        if ~fileStackTempTemp(i).isdir
            fileStackTempInd = fileStackTempInd + 1;
            strSp = strsplit(fileStackTempTemp(i).name, '_');
            fileStackTemp(fileStackTempInd).id = [strSp{1} '_' strSp{2}];
            
            fileStackTemp(fileStackTempInd).fullPathTrc = ...
                fullfile(sourceDataTrcFolder, strSp{1}, 'Only_Right_arm', [strSp{1} '_' strSp{2} '_' strSp{3} '.trc']); 
            fileStackTemp(fileStackTempInd).fullPathMat = ...
                fullfile(fileStackTempTemp(i).folder, fileStackTempTemp(i).name); 
        end
    end
    
    [cost_function_names, sortInd] = sort(cost_function_names);
    
    for i = 1:length(specStruct.ccost_array)
        specStruct.ccost_array{i} = specStruct.ccost_array{i}(sortInd);
    end
        
    runSettings.sharcNet = 0;     % is the function running on the SHARCNET cluster?
    runSettings.plotFig = 1;      % should we plot the output figures and save them?
    runSettings.saveResults = 1;  % do we need to keep a copy of the trained classifier?
    runSettings.verbose = 1;      % should we be verbose to the console or not?
    
    nowTimeStr = datestr(now, 'yyyymmdd_HHMMSS');
    outputInstancePath = [nowTimeStr];
    
%     for ind_fileStack = (0*3+3):4:length(fileStackTemp) % number of jumps        
    for ind_fileStack = 1:length(fileStackTemp)
        currInstName = fileStackTemp(ind_fileStack).id;
        fprintf('(%u/%u): %s\n', ind_fileStack, length(fileStackTemp), currInstName);

% %         if ~strcmpi(currInstName, 'Subject05_SingleArmBaseLine')
% %             fprintf('Skipping file %s\n', currInstName);
% %             continue
% %         end
        
        strsplitStr = strsplit(currInstName, '_');
        if ~strcmpi(strsplitStr(end), 'SingleArmBaseLine')
            fprintf('Skipping file %s\n', currInstName);
            continue
        end
        
        outputPath = fullfile(pathToOutput, specStruct.dataset, outputString, outputInstancePath, currInstName);
        checkMkdir(outputPath);
        
        currFilestack.dataset = specStruct.dataset;
        
        switch currFilestack.dataset                
            case {'expressive_ioc'}
                manSeg = '';
                jumpFilePath = pathToRawData;
                filesToLoad = fileStackTemp(ind_fileStack);
        end
        
        currFilestack.ccost_array = specStruct.ccost_array;
        currFilestack.cconst_array = specStruct.cconst_array;

        runSettings.nowTimeStr = nowTimeStr;
        runSettings.variableFactors = variableFactorsInput;
        
        % Find information specifying pick and place timestamps
        sequenceTime = findElementInStruct(segmentsOfInterest, currInstName);
        sequenceTime.use = 1;
         
        if calculateNormalization
            strsplitStr = strsplit(currInstName, '_');
            
            if strcmpi(strsplitStr(end), 'SingleArmBaseLine') % only run the normalization on a subset
                runSettings.variableFactors.normalizationMethod_cf = 'rang';
                param_test = setup_main(filesToLoad, [], currFilestack, run_mode, runSettings, [], [], sequenceTime);
                norm_coeff_calc = [norm_coeff_calc; param_test.coeff_cf.array];
            else
                fprintf('  Skipping...\n');
            end
            
%             average_val_use = mean(norm_coeff_calc, 1)
            
            if ind_fileStack < length(fileStackTemp) 
                continue
            else
                average_val_use = mean(norm_coeff_calc, 1);
                pause;
                continue
            end
        end
                        
        % load the data to know how much of the file needs parsing        
        [param] = setup_main(filesToLoad, '', currFilestack, 'win', runSettings, [], [], sequenceTime); % always set up in 'win' mode not 'sim'
                
        if (~isempty(fieldnames(sequenceTime)) && ~calculateNormalization) && sequenceTime.use
            startLength = round(sequenceTime.start/2);
            endLength = round(sequenceTime.end/2);
        else      
            startLength = round(length(param.t_full)/10)*10 - 30/param.dt_full + 1;
            endLength = length(param.t_full);
        end
        
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
        
        try
            windowCount = size(indToUse_window, 1);
            for ind_windowCount = 1:windowCount
                currInstNameUse = [currInstName '_' num2str(indToUse_window(ind_windowCount, 1)) '_' num2str(indToUse_window(ind_windowCount, 2))];
                variableFactors = variableFactorsInput;
                variableFactors.dataLoadLength = indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2);
                runSettings.variableFactors = variableFactors;
                
                main(filesToLoad, '', outputPath, currInstNameUse, runSettings, currFilestack, run_mode, cost_function_names, sequenceTime);
            end
            
        catch err
            errMessageTmp = regexp(err.message,',','split'); % error messages with commas in it
            errMessageFirstComma = errMessageTmp{1};
            errMessage = [errMessageFirstComma ' - ' err.stack(1).file];
            fprintf('Error START: [%s]\n', errMessage);
            
            for err_ind = 1:length(err.stack)
%                             print the error messages to screen
                err.stack(err_ind)
            end
% % %             
% % %             fprintf('Error END: [%s]\n', errMessage);
% % %         end
        end
    end
