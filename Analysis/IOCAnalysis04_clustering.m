function IOCAnalysis()
    setPaths();
    nownowstr = datestr(now, 'yyyymmddHHMMSS');
    sourceSuffix = '20200413_FatigueFull_3CF';
    targetSuffix = '20200413_FatigueFull_3CF';
    searchString = 'mat_*_3DOF_3CF*.mat';
      
    basePath = ['D:\results\fatigue_ioc03_weightsPattern\' sourceSuffix '\mat\'];
    outputPath = ['D:\results\fatigue_ioc04_weightsCluster\' targetSuffix '\'];
    
    outCsv = [outputPath 'analysis_' nownowstr '.csv'];
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
%         combinedStats = [stats_q stats_dq stats_tau stats_dtau stats_weights stats_dweights];
        combinedStats = [stats_weights stats_dweights];
        cumStats{subjectNum} = combinedStats;
    end
    
    nSubject = length(cumStats);
    allDofs = {cumStats{2}.name};
    allFeaturesSingle = cumStats{2}(1).segStats_SingleWindow.Properties.VariableNames(5:end);
    % now assemble everything into a 3D array
    for ind_subject = 1:length(cumStats)
        currSubj = cumStats{ind_subject};
        
        if isempty(currSubj)
            continue
        end
        
        for ind_dof = 1:length(allDofs)
%             currStats = cumStats{ind_subject}(ind_dof);
            currFeatureTable = cumStats{ind_subject}(ind_dof).segStats_SingleWindow;
            currIndivTable = cumStats{ind_subject}(ind_dof).regression_individual_singleWindow;
            currRegressionTable = cumStats{ind_subject}(ind_dof).regression_cumulative_singleWindow;
            
            [currFeatureTableSeg, currFeatureTableRest, segMask, restMask] = sepSegRest(currFeatureTable);
            [currRegressionTableSeg, currRegressionTableRest, segMask, restMask] = sepSegRest(currRegressionTable);
            
            for ind_features = 1:length(allFeaturesSingle)
                featureTableSegTimeSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.time;
                featureTableSegDataSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.(allFeaturesSingle{ind_features});
                
                featureTableRestTimeSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest.time;
                featureTableRestDataSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest.(allFeaturesSingle{ind_features});
            end
        end
    end
    
%     allFeaturesMultiple = cumStats{1}(1).segStats_MultipleWindow.Properties.VariableNames(8:end);
%     % now assemble everything into a 3D array
%     for ind_subject = 1:length(cumStats)
%         currSubj = cumStats{ind_subject};
%         
%         if isempty(currSubj)
%             continue
%         end
%         
%         for ind_dof = 1:length(allDofs)
%             %             currStats = cumStats{ind_subject}(ind_dof);
%             currFeatureTable = cumStats{ind_subject}(ind_dof).segStats_MultipleWindow;
%             currRegressionTable = cumStats{ind_subject}(ind_dof).segStats_MultipleWindow;
%             [currFeatureTableSeg, currFeatureTableRest, segMask, restMask] = sepSegRest(currFeatureTable);
%             [currRegressionTableSeg, currRegressionTableRest, segMask, restMask] = sepSegRest(currRegressionTable);
%             
%             for ind_features = 1:length(allFeaturesMultiple)
%                 featureTableSegTimeMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.time2;
%                 featureTableSegDataMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.(allFeaturesMultiple{ind_features});
%                 
%                 featureTableRestTimeMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableRest.time2;
%                 featureTableRestDataMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableRest.(allFeaturesMultiple{ind_features});
%             end
%         end
%     end
    
    plotStuff('SingleSeg', allDofs, allFeaturesSingle, nSubject, featureTableSegTimeSingle, featureTableSegDataSingle, outCsv, outputPath);
    plotStuff('SingleRest', allDofs, allFeaturesSingle, nSubject, featureTableRestTimeSingle, featureTableRestDataSingle, outCsv, outputPath);
%     plotStuff('MultipleSeg', allDofs, allFeaturesMultiple, nSubject, featureTableSegTimeMultiple, featureTableSegDataMultiple, outCsv, outputPath);
%     plotStuff('MultipleRest', allDofs, allFeaturesMultiple, nSubject, featureTableRestTimeMultiple, featureTableRestDataMultiple, outCsv, outputPath);
end

function plotStuff(typeLabel, allDofs, allFeaturesSingle, nSubject, featureTableSegTime, featureTableSegData, outCsv, outputPath) 
    subjectColours = distinguishable_colors(nSubject);
    currInd = 0;
    
    % then plot everything
    for ind_dof = 1:length(allDofs)
        for ind_features = 1:length(allFeaturesSingle)
            currInd = currInd + 1;
            fprintf('Currently on %s - %u of %u\n', typeLabel, currInd, length(allDofs)*length(allFeaturesSingle));
            
            currDof = allDofs{ind_dof};
            currFeature = allFeaturesSingle{ind_features};
            figName = ['combined_' typeLabel '_' currDof '_' currFeature];
            figSavePath = fullfile(outputPath, figName);
            
            figName2 = ['indiv_' typeLabel '_' currDof '_' currFeature];
            figSavePath2 = fullfile(outputPath, figName2);
            
            h1 = figure('Position', [-1919 69 1920 964.8000]);
            h2 = figure('Position', [-1919 69 1920 964.8000]);
            hold on; 
            
            bSign = zeros(3, 1);
            
            for ind_subjects = 1:nSubject
                currColour = subjectColours(ind_subjects, :);
                currTime = featureTableSegTime{ind_dof, ind_features, ind_subjects};
                currData = featureTableSegData{ind_dof, ind_features, ind_subjects};
                
                if isempty(currTime)
                    continue
                end
                
                [b(:, ind_subjects), Rsq2(ind_subjects), X, yCalc2] = linearFit(currTime, currData);
                dataLen(ind_subjects) = length(currData);
                
                signageB = sign(b(2, ind_subjects));
                if signageB > 0
                    bSign(1) = bSign(1) + 1;
                elseif signageB < 0
                    bSign(2) = bSign(2) + 1;
                elseif isnan(b(2, ind_subjects))
                    bSign(3) = bSign(3) + 1;
                end
                
                currLabel = ['S' num2str(ind_subjects) ', b=' num2str(b(2, ind_subjects), '%0.4f'), ', R2=', num2str(Rsq2(ind_subjects), '%0.4f')];
                
                figure(h1);
                hold on;
                plot(currTime, currData, 'DisplayName', currLabel, 'Color', currColour, 'MarkerSize', 14, 'Marker', 'o', 'LineStyle', 'none');
                ph = plot(currTime, yCalc2, 'Color', currColour);
                ph.Annotation.LegendInformation.IconDisplayStyle = 'off';
                
                figure(h2);
                subplot(3, 5, ind_subjects);
                hold on;
                plot(currTime, currData, 'DisplayName', currLabel, 'Color', currColour, 'MarkerSize', 14, 'Marker', 'o', 'LineStyle', 'none');
                ph = plot(currTime, yCalc2, 'Color', [1 0 0]);
                ph.Annotation.LegendInformation.IconDisplayStyle = 'off';
                title(currLabel);
            end
            
            % what is the mean/std slope and Rsq value? 
            meanB = mean(b(2, :));
            stdB = std(b(2, :));
            meanRsq = mean(Rsq2);
            stdRsq = std(Rsq2);
            
            figure(h1);
            titleStr = [typeLabel '_' currDof '_' currFeature, ...
                ', b=' num2str(meanB, '%0.4f') '\pm' num2str(stdB, '%0.4f'), ...
                ', bsign=(+)' num2str(bSign(1), '%0.4f') ', (-)' num2str(bSign(2), '%0.4f'), ...
                ', R2=', num2str(meanRsq, '%0.4f') '\pm' num2str(stdRsq, '%0.2f')];
            title(titleStr);
            
            legend('show');
            
            saveas(h1, figSavePath, 'png');
            saveas(h1, figSavePath, 'fig');
            close(h1);
            
            saveas(h2, figSavePath2, 'png');
            saveas(h2, figSavePath2, 'fig');
            close(h2);
            
            if ~exist(outCsv, 'file')
                header = 'typeLabel,dof,feature,b_mean,b_std,R2_mean,R2_std,bsign_plus,bsign_minus,bsign_nan';
                for ind_subjects = 1:nSubject
                    header = [header ',bmag_s' num2str(ind_subjects)];
                end
                for ind_subjects = 1:nSubject
                    header = [header ',r2_s' num2str(ind_subjects)];
                end
                for ind_subjects = 1:nSubject
                    header = [header ',len_s' num2str(ind_subjects)];
                end
                header = [header '\n'];
            else
                header = '';
            end
            
            fid = fopen(outCsv, 'a');
            fprintf(fid, [header]);
            
            fprintf(fid, '%s,%s,%s', typeLabel, currDof, currFeature);
            fprintf(fid, ',%f,%f', meanB, stdB);
            fprintf(fid, ',%f,%f', meanRsq, stdRsq);
            fprintf(fid, ',%f,%f,%f', bSign(1), bSign(2), bSign(3));
            
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', b(2, ind_subjects));
            end
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', Rsq2(ind_subjects));
            end
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', dataLen(ind_subjects));
            end
            fprintf(fid, '\n');
            fclose(fid);
        end
    end    
end

