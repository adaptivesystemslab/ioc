function IOCAnalysis()
    setPaths();
    
%     perturbAmount = [1e-1 1e-2 1e-3 1e-4];
%     perturbAmount = [1e-1 1e-2 1e-3 1e-4]*1e4;
%     residualThreshold = [1e1 1e0 1e-1 1e-2 1e-3 1e-4 1e-5]; % 0.1
%     residualThreshold = [1e1 1e0 1e-1 1e-2 1e-3 1e-4 1e-5]/1e7; 

%     perturbAmount = 10.^[-5:4]; 
%     residualThreshold = 10.^[-7:10];
    perturbAmount = [1e-1];
    residualThreshold = [1e-1];

    maxAcross = 4;
    maxDown = 7;
    
    faceColours = brewermap(8, 'paired');
    
    nowstr = datestr(now, 'yyyymmddHHMMSS');
    inputBasePath = 'D:\aslab_github\ioc\Data\IOC\';
    outputBasePath = ['D:\aslab_github\ioc\Data\IOC\plot_' nowstr];
    masterPathCsv = [outputBasePath '\a_summary.csv'];
    checkMkdir(outputBasePath);
    
    % iterate through all the folders in bathpath to produce figures and
    % table summary of the results
    
    % look for 'weights*1_' since that will be unique
    % then
    
    dirPath = dir(inputBasePath);
    
    for i = 1:length(dirPath)
        currPath = fullfile(inputBasePath, dirPath(i).name);
        
        if strcmpi(currPath(end), '.')
            continue; % it's . or ..
        end
        
        dirCurrSubPath = dir(currPath);
        
        for j = 1:length(dirCurrSubPath)
            currSubPath = fullfile(currPath, dirCurrSubPath(j).name);
            if strcmpi(currSubPath(end), '.')
                continue; % it's . or ..
            end
            
            if ~exist(currSubPath, 'dir')
                continue;
            end
            
%             if ~strcmpi(dirCurrSubPath(j).name, 'Subj1_7DOF_8CF')
%                 continue;
%             end
            
            matData = assembleData(currSubPath);
            
            if isempty(matData)
                continue;
            end
            
            suffix = [dirPath(i).name '_' matData.trialInfo.runName '_' matData.trialInfo.templateName];
            
            outputPath = [outputBasePath '\' dirPath(i).name];
            checkMkdir(outputPath);
%             try
                % plot results
%                 fprintf('%s\n', suffix);
                outputPathFig1 = fullfile(inputBasePath, ['fig_results_individual_' suffix]);
                outputPathFig2 = fullfile(outputPath, ['fig_results_cumulativeAllPass_' suffix]);
                outputPathFig3 = fullfile(outputPath, ['fig_results_cumulativeRankPass_' suffix]);
                outputPathMat1 = fullfile(outputPath, ['mat_results_cumulativeAllPass_' suffix]);
                outputPathCsv = fullfile(outputPath, ['csv_' suffix]);
%                 csv_populate(matData, masterPathCsv);
                plotting_individual(matData, outputPathFig1, outputPathCsv, masterPathCsv);
                matSave = plotting_cumulative(matData, outputPathFig2, outputPathFig3, outputPathCsv, masterPathCsv, faceColours, outputPathMat1);
                
% % %                 % load data and perturb
% % %                 for ind_perturbAmount = 1:length(perturbAmount)
% % %                     for ind_residual = 1:length(residualThreshold)
% % %                         suffix_suffix = [num2str(ind_perturbAmount) '_' num2str(ind_residual)];
% % %                         suffix2 = [suffix '_' suffix_suffix];
% % %                         outputPathFig1 = fullfile(outputPath, ['fig_perturb_single_' suffix2]);
% % %                         outputPathFig2 = fullfile(outputPath, ['fig_perturb_combined_' suffix2]);
% % %                         suffix3 = [suffix '_' num2str(ind_perturbAmount) '_' num2str(ind_residual)]
% % %                         outputPathFig1 = fullfile(outputPath, ['fig_perturb_single_' suffix3]);
% % %                         outputPathFig2 = fullfile(outputPath, ['fig_perturb_combined_' suffix3]);
% % %                         matFilePath = fullfile(outputPath, ['mat_perturb_combined_' suffix3]);
% % %                         %                        peturbWeights_single(matData, outputPathFig1, currPerturb, currResidual);
% % %                         
% % %                         fprintf('Processing %s, %s, %s\n', suffix, suffix2, suffix3);
% % %                         
% % %                         [t{ind_perturbAmount, ind_residual}, weights_cum_cum{ind_perturbAmount, ind_residual}, weights_blocks{ind_perturbAmount, ind_residual}, ...
% % %                             weight_labels{ind_perturbAmount, ind_residual}, residual_keep_cumulative{ind_perturbAmount, ind_residual}] = ...
% % %                             peturbWeights_combined(matData,d outputPathFig2, perturbAmount(ind_perturbAmount), residualThreshold(ind_residual), faceColours, suffix_suffix, ...
% % %                             matSave, matFilePath);
% % %                     end
% % %                 end
% % %                 
% % %                 figFileName = fullfile(outputPath, ['a_fig_perturb_combinedGridWeight_' suffix]);
% % %                 plotWeights_Grid(figFileName, perturbAmount, residualThreshold, maxAcross, maxDown, faceColours,  ...
% % %                     t, weights_cum_cum, residual_keep_cumulative, [], weight_labels);
% % %                 
% % %                 figFileName = fullfile(outputPath, ['a_fig_perturb_combinedGridRemoval_' suffix]);
% % %                 plotWeights_Grid(figFileName, perturbAmount, residualThreshold, maxAcross, maxDown, faceColours,  ...
% % %                     t, [], residual_keep_cumulative, weights_blocks, weight_labels);
                       
%             catch err
%                 err
%                 close all;
%             end
        end
    end
end