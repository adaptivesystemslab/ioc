function IOCAnalysis()
    setPaths();
%     nowstr = datestr(now, 'yyyymmddHHMMSS');
    nowstr = '20200316_fatigueEdges';
    nowstr2 = '20200316_fatigueEdges';
      
    basePath = ['D:\results\fatigue_ioc03_weightsPattern\' nowstr '\mat\'];
    searchString = 'mat_*_3DOF_3CF*.mat';
    outputPath = ['D:\results\fatigue_ioc04_weightsCluster\' nowstr2 '\'];
    outCsv = [outputPath 'analysis.csv'];
    checkMkdir(outputPath);
    
    currBasePathDir = dir([basePath searchString]);
    for j = 1:length(currBasePathDir)
        currFileName = currBasePathDir(j).name;
        currFullPath = fullfile(basePath, currFileName);
        
        if strcmpi(currFileName(end), '.')
            continue; % it's . or ..
        end

        % load each entry and prepare to plot
        load(currFullPath);
        
        subjectNum = str2num(trialInfo.runName(8:9));
        combinedStats = [stats_q stats_dq stats_tau stats_weights];
        cumStats{subjectNum} = combinedStats;
    end
    
    nSubject = length(cumStats);
    allDofs = {cumStats{1}.name};
    allFeatures = cumStats{1}(1).segStats_SingleWindow.Properties.VariableNames(5:end);
    
    % now assemble everything into a 3D array
    for ind_subject = 1:length(cumStats)
        currSubj = cumStats{ind_subject};
        
        if isempty(currSubj)
            continue
        end
        
        for ind_dof = 1:length(allDofs)
%             currStats = cumStats{ind_subject}(ind_dof);
            currFeatureTable = cumStats{ind_subject}(ind_dof).segStats_SingleWindow;
            currRegressionTable = cumStats{ind_subject}(ind_dof).regression_singleWindow;
            [currFeatureTableSeg, currFeatureTableRest, segMask, restMask] = sepSegRest(currFeatureTable);
            [currRegressionTableSeg, currRegressionTableRest, segMask, restMask] = sepSegRest(currRegressionTable);
            
            for ind_features = 1:length(allFeatures)
                featureTableSegTime{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.time;
                featureTableSegData{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.(allFeatures{ind_features});
                
                featureTableRestTime{ind_dof, ind_features, ind_subject} = currFeatureTableRest.time;
                featureTableRestData{ind_dof, ind_features, ind_subject} = currFeatureTableRest.(allFeatures{ind_features});
            end
        end
    end
    
    subjectColours = distinguishable_colors(nSubject);
    
    % then plot everything
    for ind_dof = 1:length(allDofs)
        for ind_features = 1:length(allFeatures)
            currDof = allDofs{ind_dof};
            currFeature = allFeatures{ind_features};
            figName = [currDof '_' currFeature '_seg'];
            figSavePath = fullfile(outputPath, figName);
            
            h = figure('Position', [-1919 69 1920 964.8000]);
            hold on; 
            
            bSign = zeros(2, 1);
            
            for ind_subjects = 1:nSubject
                currColour = subjectColours(ind_subjects, :);
                currTime = featureTableSegTime{ind_dof, ind_features, ind_subjects};
                currData = featureTableSegData{ind_dof, ind_features, ind_subjects};
                
                if isempty(currTime)
                    continue
                end
                
                [b(:, ind_subjects), Rsq2(ind_subjects), X, yCalc2] = linearFit(currTime, currData);
                
                signageB = sign(b(2, ind_subjects))
                if signageB > 0
                    bSign(1) = bSign(1) + 1;
                elseif signageB < 0
                    bSign(2) = bSign(2) + 1;
                end
                
                
                currLabel = ['S' num2str(ind_subjects) ', b=' num2str(b(2, ind_subjects), '%0.4f'), ', R2=', num2str(Rsq2(ind_subjects), '%0.4f')];
                
                plot(currTime, currData, 'DisplayName', currLabel, 'Color', currColour, 'MarkerSize', 14, 'Marker', 'o', 'LineStyle', 'none');
                ph = plot(currTime, yCalc2, 'Color', currColour);
                ph.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            
            % what is the mean/std slope and Rsq value? 
            meanB = mean(b(2, :));
            stdB = std(b(2, :));
            meanRsq = mean(Rsq2);
            stdRsq = std(Rsq2);
            
            titleStr = [currDof '_' currFeature, ...
                ', b=' num2str(meanB, '%0.4f') '\pm' num2str(stdB, '%0.4f'), ...
                ', bsign=(+)' num2str(bSign(1), '%0.4f') ', (-)' num2str(bSign(2), '%0.4f'), ...
                ', R2=', num2str(meanRsq, '%0.4f') '\pm' num2str(stdRsq, '%0.2f')];
            title(titleStr);
            
            legend('show');
            
            saveas(h, figSavePath, 'png');
            saveas(h, figSavePath, 'fig');
            close(h);
            
            if ~exist(outCsv, 'file')
                header = 'dof,feature,b_mean,b_std,bsign_plus,bsign_minus,R2_mean,R2_std';
                header = [header '\n'];
            else
                header = '';
            end
            
            fid = fopen(outCsv, 'a');
            fprintf(fid, [header]);
            
            fprintf(fid, '%s,%s', currDof, currFeature);
            fprintf(fid, ',%f,%f', meanB, stdB);
            fprintf(fid, ',%f,%f', bSign(1), bSign(2));
            fprintf(fid, ',%f,%f', meanRsq, stdRsq);
            
            fprintf(fid, '\n');
            fclose(fid);
        end
    end
 
    
%      h = figure('Position', [-1919 69 1920 964.8000]);
%         hold on
    
%     for ind_dof = 1:length(cumStats{1})
%        
%         
%         data = [];
%         for ind_subject = 1:length(cumStats)
%             c
% 
%             for ind_segment = 1:length(currStats.segmentStatsSeg)
%                 data = [data currStats.segmentStatsSeg{ind_segment}.];
%             end
%             
%         end
%     end
end

