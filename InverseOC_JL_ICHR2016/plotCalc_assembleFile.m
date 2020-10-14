 % assemble the final string over the various files
    indToUse_window = [];
    output_inverse = {};
    feature_win_save = {};
    feature_recon_local = {};
    t_recon_plot_array = {};
    q_recon_plot_array = {};
    dq_recon_plot_array = {};
    ddq_recon_plot_array = {};
    dqtau_recon_plot_array = {};
    ddx_recon_plot_array = {};
    dqtau_plot_array = {};
    
    minRmseIndArray = [];
    segmentInfoCount = [];
    segmentInfoStart = [];
    segmentInfoEnd = [];
    elapsedTimeTotal = 0;
    minRmseInd_plot = [];
    dofsUsedInIOC_plot = [];
    for ind = 1:length(matFilesToLoad)
        % load each individual file
        fprintf('  Loading window %u of %u\n', ind, length(matFilesToLoad));
        
        currFolder = pwd;
        [pathstr,name,ext] = fileparts(matFilesToLoad{ind});
        cd(pathstr);
        currLoadPackage_tmp = load([name ext]);
		currLoadPackage = currLoadPackage_tmp.saveVar;
        cd(currFolder);
        
        if ind < length(matFilesToLoad)
            currIndRange = indTotal(ind, :);
            currIndRange(2) = indTotal(ind+1, 1);
        else
            currIndRange = indTotal(ind, :);
        end
        
        % pull out the window range
        if currLoadPackage.indToUse_window(1, 1) == 1 
            indToUse_windowToAdd = currLoadPackage.indToUse_window(:, 1) + currIndRange(1) - 1;
            startWindowIndToAdd = find(indToUse_windowToAdd == currIndRange(1));
            endWindowIndToAdd = find(indToUse_windowToAdd == currIndRange(2)) - 1; % -1 to prevent overlap with next loaded window
        else
            indToUse_windowToAdd = currLoadPackage.indToUse_window(:, 1) + currIndRange(1) - 1;
            startWindowIndToAdd = 1;
            endWindowIndToAdd = size(currLoadPackage.indToUse_window, 1);
        end
        
        if isempty(endWindowIndToAdd) || endWindowIndToAdd == 0
            endWindowIndToAdd = length(indToUse_windowToAdd);
        end
        
        targetWindowToPull = startWindowIndToAdd:endWindowIndToAdd;
        
        targetWindowIndToHit = currIndRange(1):currIndRange(2);
        targetWindowIndToPull = [currIndRange(1):currIndRange(2)] - currIndRange(1) + 1;
        feature_full.t(:, targetWindowIndToHit) = currLoadPackage.feature_full.t(:, targetWindowIndToPull);
        feature_full.q(:, targetWindowIndToHit) = currLoadPackage.feature_full.q(:, targetWindowIndToPull);
        feature_full.dq(:, targetWindowIndToHit) = currLoadPackage.feature_full.dq(:, targetWindowIndToPull);
        feature_full.ddq(:, targetWindowIndToHit) = currLoadPackage.feature_full.ddq(:, targetWindowIndToPull);
        ccost_full(:, targetWindowIndToHit) = currLoadPackage.ccost_full(targetWindowIndToPull, :)';
        
        packageToAdd = currLoadPackage.output_inverse(targetWindowToPull);
        output_inverse = [output_inverse packageToAdd];
        
        packageToAdd = currLoadPackage.t_recon_plot_array(targetWindowToPull);
        t_recon_plot_array = [t_recon_plot_array packageToAdd];
        
        packageToAdd = currLoadPackage.elapsedTime;
        elapsedTimeTotal = elapsedTimeTotal + packageToAdd;
                
        packageToAdd = currLoadPackage.q_recon_plot_array(targetWindowToPull);
        q_recon_plot_array = [q_recon_plot_array packageToAdd];
        
        packageToAdd = currLoadPackage.dq_recon_plot_array(targetWindowToPull);
        dq_recon_plot_array = [dq_recon_plot_array packageToAdd];
        
        packageToAdd = currLoadPackage.ddq_recon_plot_array(targetWindowToPull);
        ddq_recon_plot_array = [ddq_recon_plot_array packageToAdd];
        
%         for j = 1:length(targetWindowToPull)
%             currMin = currLoadPackage.minRmseInd_plot(targetWindowToPull(j));
% %             packageToAdd = currLoadPackage.feature_recon_local{targetWindowToPull(j)}{currMin}.q;
% %             q2_recon_plot_array = [q2_recon_plot_array packageToAdd];
%             
%             packageToAdd = currLoadPackage.feature_recon_local{targetWindowToPull(j)}{currMin}.en;
%             dqtau_recon_plot_array = [dqtau_recon_plot_array packageToAdd];
%             
%             packageToAdd = currLoadPackage.feature_recon_local{targetWindowToPull(j)}{currMin}.ddx;
%             ddx_recon_plot_array = [ddx_recon_plot_array packageToAdd];
%             
%             packageToAdd = currLoadPackage.feature_win_save{targetWindowToPull(j)}.en;
%             dqtau_plot_array = [dqtau_plot_array packageToAdd];
%         end
        
        packageToAdd = currLoadPackage.minRmseIndArray(targetWindowToPull);
        minRmseIndArray = [minRmseIndArray packageToAdd];
        
        packageToAdd = currLoadPackage.feature_win_save(targetWindowToPull);
        feature_win_save = [feature_win_save packageToAdd];
        
        packageToAdd = currLoadPackage.feature_recon_local(targetWindowToPull);
        feature_recon_local = [feature_recon_local packageToAdd];
        
        packageToAdd = currLoadPackage.indToUse_window(targetWindowToPull, :) + currIndRange(1) - 1;
        indToUse_window = [indToUse_window; packageToAdd];
        
        if strcmpi(currLoadPackage.run_mode, 'sim')
            packageToAdd = 1:length(currLoadPackage.segmentInfo.timeStart);
        else
            packageToAdd = 1:length(currLoadPackage.segmentInfo.timeStart);
        end
        segmentInfoCount = [segmentInfoCount; packageToAdd];
        
        packageToAdd = currLoadPackage.segmentInfo.timeStart;
        segmentInfoStart = [segmentInfoStart; packageToAdd];
        
        packageToAdd = currLoadPackage.segmentInfo.timeEnd;
        segmentInfoEnd = [segmentInfoEnd; packageToAdd];
        
        packageToAdd = currLoadPackage.minRmseInd_plot(targetWindowToPull);
        minRmseInd_plot = [minRmseInd_plot; packageToAdd];
    end
    
    [segmentCountUnique, segmentCountInd] = unique(segmentInfoCount);
    segmentInfo.segmentCount = segmentInfoCount(segmentCountInd);
    segmentInfo.timeStart = segmentInfoStart(segmentCountInd);
    segmentInfo.timeEnd = segmentInfoEnd(segmentCountInd);
    
    param = currLoadPackage.param;
    cost_function_names = currLoadPackage.cost_function_names;
    runSettings = currLoadPackage.runSettings;
%     lenToUseT = currLoadPackage.lenToUseT;
    run_mode = currLoadPackage.run_mode;
    currFilestack = currLoadPackage.currFilestack;
    
%     svd_ratio_threshold = runSettings.correlation_threshold;
    
    plotType = 'weight';
    
    currDataset = currLoadPackage.currFilestack.dataset;
    clear currLoadPackage;