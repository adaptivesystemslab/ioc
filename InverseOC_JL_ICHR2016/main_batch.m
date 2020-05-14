function main_batch(specStruct, run_mode, cost_function_names, outputString, variableFactorsInput, fileStackTemp, pathToOutput)
    % preamble loading of background parameters
    
    addpath(genpath(fullfile('Common')));
    addpath(genpath('Libraries/Robotics_Corke'));
    addpath(genpath('Libraries/rl/ik_framework/common'));
    addpath(genpath('Libraries/rl/ik_framework/instance_expressive'));
    addpath(genpath('Libraries/rl/General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));
    
    calculateNormalization = 0;
    
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

        outputPath = fullfile(pathToOutput, specStruct.dataset, outputString, outputInstancePath, currInstName);
        checkMkdir(outputPath);
        
        currFilestack.dataset = specStruct.dataset;
        
        switch currFilestack.dataset                
            case {'expressive_ioc'}
                filesToLoad = fileStackTemp(ind_fileStack);
        end
        
        currFilestack.ccost_array = specStruct.ccost_array;
        currFilestack.cconst_array = specStruct.cconst_array;

        runSettings.nowTimeStr = nowTimeStr;
        runSettings.variableFactors = variableFactorsInput;
        
        % Find information specifying pick and place timestamps
         
        if calculateNormalization
            strsplitStr = strsplit(currInstName, '_');
            
            if strcmpi(strsplitStr(end), 'SingleArmBaseLine') % only run the normalization on a subset
                runSettings.variableFactors.normalizationMethod_cf = 'rang';
                param_test = setup_main(filesToLoad, [], currFilestack, run_mode, runSettings, [], [], []);
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
        [param] = setup_main(filesToLoad, '', currFilestack, 'win', runSettings, [], [], []); % always set up in 'win' mode not 'sim'
                 
        startLength = 1;
        endLength = 5000;
%         endLength = round(length(param.t_full)/10)*10; 
%         endLength = length(param.t_full);
        
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
                
                main(filesToLoad, '', outputPath, currInstNameUse, runSettings, currFilestack, run_mode, cost_function_names, []);
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
