function [measRmseId1, measRmseId2, sensorsUsed, sensorsTotal] = plot_measurements(...
    featureSet_orig, featureSet_post, measurement_orig, measurement_post, str_orig, str_post, filepathFig, ...
    ekf_markerMask, ekf_eventhandler, ekf_markerTally)

    rmseFct = setRmseFct();
    neFct = setNormErrorFct();

    colours = distinguishable_colors(6);
    measurement_labels = featureSet_orig.measurement_labels;
    lenJointLabels = length(measurement_labels);

    measRmseId1 = zeros(1, lenJointLabels);
    measRmseId2 = zeros(1, lenJointLabels);
    sensorsTotal = zeros(1, lenJointLabels);
    sensorsUsed = zeros(1, lenJointLabels);

    for i = 1:lenJointLabels
        orig_time = featureSet_orig.time - featureSet_orig.time(1);
        orig_mes = measurement_orig(:, i).getMesArray;
        post_time = featureSet_post.time - featureSet_post.time(1);
        post_mes = measurement_post(:, i).getMesArray;
        if ~isempty(ekf_markerMask)
            orig_mask = ekf_markerMask.swappedmissing(:, i);
        else
            orig_mask = ones(size(orig_time));
        end
        
        switch size(post_mes, 2)
            case 1
                indsId1 = 1;
                indsId2 = [];
                
            case 3
                indsId1 = 1:3;
                indsId2 = [];
                
            case 6
                indsId1 = 1:3;
                indsId2 = 4:6;
        end
        
        measRmseId1(i) = rmseFct(orig_mes(orig_mask == 1, indsId1), post_mes(orig_mask == 1, indsId1));
        rmseTraj1(:, i) = neFct(orig_mes(:, indsId1), post_mes(:, indsId1));
        
        if isempty(indsId2)
            measRmseId2(i) = rmseFct(orig_mes(orig_mask == 1, indsId2), post_mes(orig_mask == 1, indsId2));
            rmseTraj2(:, i) = neFct(orig_mes(:, indsId2), post_mes(:, indsId2));
        end
    end

    for i = 1:lenJointLabels
        orig_time = featureSet_orig.time - featureSet_orig.time(1);
        orig_mes = measurement_orig(:, i).getMesArray;
        post_time = featureSet_post.time - featureSet_post.time(1);
        post_mes = measurement_post(:, i).getMesArray;
%         orig_mask = featureSet_orig.measurement_mask(:, i);

        if ~isempty(ekf_markerMask)
            mask_swap = ekf_markerMask.swapped(:, i);
            mask_missing = ekf_markerMask.missing(:, i);
            
            usePercent = ekf_markerTally.correct(i)/length(orig_time);
            missingPercent = ekf_markerTally.missing(i)/length(orig_time);
            swapPercent = ekf_markerTally.swapped(i)/length(orig_time);
        else
            mask_swap = ones(size(orig_time));
            mask_missing = ones(size(orig_time));
            
            usePercent = 1;
            missingPercent = 0;
            swapPercent = 0;
        end
        
        lenDofs = size(post_mes, 2);
        
        h = figure('position', 1.0e+03 *[ 0.0010    0.0410    1.5360    0.7488]);

        useStr = [' (use=' num2str(usePercent) ',missing=' num2str(missingPercent) ',swap=' num2str(swapPercent) ')'];
        
        if size(post_mes, 2) >= 6
            ax(1) = subplot(211);
            for j = 1:3
                plot(orig_time, orig_mes(:, j), 'DisplayName', str_orig, 'Color', colours(j, :), 'LineStyle', '--');
                hold on;
                plot(post_time, post_mes(:, j), 'DisplayName', str_post, 'Color', colours(j, :), 'LineStyle', '-');
            end
%             plot(orig_time, rmseTraj(:, 1), 'DisplayName', 'norm(x-y)', 'Color', [0 0 0], 'LineStyle', '-.');
%             plot(orig_time(mask_swap == 0), rmseTraj(mask_swap == 0, 1), 'Color', [0 0 1], 'DisplayName', 'Unused markers', 'LineStyle', 'none', 'Marker', 'x');
%             plot(orig_time(mask_missing == 0), rmseTraj(mask_missing == 0, 1), 'Color', [0 1 0], 'DisplayName', 'Unused markers', 'LineStyle', 'none', 'Marker', 'x');
     
            titleStr = [measurement_labels{i} '_sensorType' num2str(measurement_orig(1, i).type) '_ID1' ...
                ', rmse=' num2str(measRmseId1(i) ) ')' useStr];
            title(titleStr);

            ax(2) = subplot(212);
            for j = 4:6
                plot(orig_time, orig_mes(:, j), 'DisplayName', str_orig, 'Color', colours(j, :), 'LineStyle', '--');
                hold on;
                plot(post_time, post_mes(:, j), 'DisplayName', str_post, 'Color', colours(j, :), 'LineStyle', '-');
            end
%             plot(orig_time, rmseTraj(:, 2), 'DisplayName', 'norm(x-y)', 'Color', [0 0 0], 'LineStyle', '-.');
%             plot(orig_time(orig_mask == 0), rmseTraj(orig_mask == 0, 2), 'Color', [0 0 0], 'DisplayName', 'Unused markers', 'Marker', 'x');

            titleStr = [measurement_labels{i} '_sensorType' num2str(measurement_orig(1, i).type) '_ID2' ...
                ', rmse=' num2str(measRmseId2(i) ) ')' useStr];
            title(titleStr);
            
            linkaxes(ax, 'xy');
        else
            for j = 1:lenDofs
                hold on;
                plot(orig_time, orig_mes(:, j), 'DisplayName', str_orig, 'Color', colours(j, :), 'LineStyle', '-');
                plot(post_time, post_mes(:, j), 'DisplayName', str_post, 'Color', colours(j, :), 'LineStyle', '--');
            end
            
            plot(orig_time, rmseTraj1(:, i), 'DisplayName', 'norm(x-y)', 'Color', [0 0 0], 'LineStyle', '-.');
            plot(orig_time(mask_swap == 0), zeros(size(orig_time(mask_swap == 0))), 'Color', [0.3 0 0.3], 'DisplayName', 'Swapped markers', 'LineStyle', 'none', 'Marker', 'o');
            plot(orig_time(mask_missing == 0), zeros(size(orig_time(mask_missing == 0))), 'Color', [0.3 0.3 0], 'DisplayName', 'Missing markers', 'LineStyle', 'none', 'Marker', 'x');
    
            titleStr = [measurement_labels{i} '_sensorType' num2str(measurement_orig(1, i).type) '_ID1' ...
                ', rmse=' num2str(measRmseId1(i) ) ' ' useStr];
            title(titleStr);
        end
        
        legend show;
        excludeVal = 1e5-5;
        ymin = floor(min(min([orig_mes(orig_mes < excludeVal); post_mes(post_mes < excludeVal)])));
        ymax = ceil(max(max([orig_mes(orig_mes < excludeVal); post_mes(post_mes < excludeVal)])));
        if ymin ~= ymax
            ylim([ymin ymax]);
        end
        
        finalFilePath = [filepathFig '_' measurement_labels{i}];
        saveas(h, finalFilePath, 'png');
        saveas(h, finalFilePath, 'fig');

        close(h);
    end
    
%     % hardwired marker data
%     markersOfInterest = {...
%         {'ASIS_R', 'KNEE_R_IMU_T'}, ...
%         {'KNEE_R_MED', 'KNEE_R_IMU_T'}, ...
%         {'KNEE_R_MED', 'ANKLE_R_IMU_T'}, ...
%         {'ANKLE_R_MED', 'ANKLE_R_IMU_T'}};
%     
%     for i = 1:length(markersOfInterest)
%         orig_time = featureSet_orig.time - featureSet_orig.time(1);
%          
%         indSource = find(ismember(featureSet_orig.measurement_labels, markersOfInterest{i}{1}));
%         indTarget = find(ismember(featureSet_orig.measurement_labels, markersOfInterest{i}{2}));
%         orig_source = measurement_orig(:, indSource).getMesArray;
%         orig_target = measurement_orig(:, indTarget).getMesArray;
%         orig_mes = orig_source - orig_target;
%         post_source = measurement_post(:, indSource).getMesArray;
%         post_target = measurement_post(:, indTarget).getMesArray;
%         post_mes = post_source - post_target;
%         
%         lenDofs = size(post_mes, 2);
%         
%          h = figure('position', 1.0e+03 *[ 0.0010    0.0410    1.5360    0.7488]);
%         for j = 1:lenDofs
%             hold on;
%             plot(orig_time, orig_mes(:, j), 'DisplayName', str_orig, 'Color', colours(j, :), 'LineStyle', '-');
%             plot(post_time, post_mes(:, j), 'DisplayName', str_post, 'Color', colours(j, :), 'LineStyle', '--');
%         end
%         
%         legend show;
%         excludeVal = 1e5-5;
%         ymin = floor(min(min([orig_mes(orig_mes < excludeVal); post_mes(post_mes < excludeVal)])));
%         ymax = ceil(max(max([orig_mes(orig_mes < excludeVal); post_mes(post_mes < excludeVal)])));
%         if ymin ~= ymax
%             ylim([ymin ymax]);
%         end
%         
%         titleStr = [markersOfInterest{i}{1} '_' markersOfInterest{i}{2}];
%         title(titleStr);
%         
%         finalFilePath = [filepathFig '_' markersOfInterest{i}{1} '_' markersOfInterest{i}{2} '_offset'];
%         saveas(h, finalFilePath, 'png');
%         saveas(h, finalFilePath, 'fig');
% 
%         close(h);
%     end
end