function batchConvertClustersToSegmentPoints
    %     sysparam.tolerance = 0.2;
    %     sysparam.offset = 0;
    %     sysparam.offsetDir = 'shift';

%     offsetArray1 = [0];
%     offsetArray2 = [-0.05];
        offsetArray1 = [-0.06:0.01:0.01];
        offsetArray2 = [-0.06:0.01:0.01];

    % generate paths to the .mat files
    sourcePathRoot = 'D:\MATLABResults\STAT841Proj-Healthy1\script20_healthy1-All_4_SVMradialonly1\';
    sourceSuffix = '15_4_110_PCA_0_SVM_radial_None';
    exportPath = fullfile(sourcePathRoot, ['clusterdata_' sourceSuffix '.csv']);
    
    fileStack = {};
    fileStackCounter = 0;
    sourceRootDir = dir(sourcePathRoot);
    for i = 1:length(sourceRootDir)
        if ~sourceRootDir(i).isdir
            continue
        elseif length(sourceRootDir(i).name) < 15
            continue
        end
        
        sourceSubPath = fullfile(sourcePathRoot, sourceRootDir(i).name, sourceSuffix, 'Testing');
        sourceSubDir = dir(fullfile(sourceSubPath, '*.mat'));
        for i2 = 1:length(sourceSubDir)
            fileStackCounter = fileStackCounter + 1;
            fileStack{fileStackCounter} = fullfile(sourceSubPath, sourceSubDir(i2).name);
        end
    end

    % load all the files in this folder

    arrayLength = length(offsetArray1)*length(offsetArray2);
    maxIterLength = length(offsetArray1)*length(offsetArray2)*length(fileStack);
    summaryMatrix{arrayLength} = []; 
    assessConfig(arrayLength, 3) = 0;
    
    progressBarCounter = 0;
    for i = 1:length(fileStack)
%         sourcePath = 'D:\MATLABResults\STAT841Proj-Healthy1\2014-03-07-12-36-50\2014-03-07-train_5-test_1\15_25_110_PCA_0_SVM_radial_None\Subject1_KEFO_SIT_SLO1.mat';
        
        % load source data
        load(fileStack{i});    
        [manualSegStruct, algSegStruct] = convertClustersToSegmentPoints(saveStruct);
%         close all
%         continue
        counter = 0;
        for ind_offsetArray1 = 1:length(offsetArray1)
            fprintf('Currently in %u/%u \n', progressBarCounter, maxIterLength);
            
            for ind_offsetArray2 = 1:length(offsetArray2)
                progressBarCounter = progressBarCounter + 1;
                    
                sysparam.offsetDir = 'expand'; % 'shift';
                sysparam.offset = [offsetArray1(ind_offsetArray1) offsetArray2(ind_offsetArray2)];
                sysparam.tolerance = 0.2;
                
                metric = manualSegmentationComparison(manualSegStruct, algSegStruct, sysparam);
                
                counter = counter + 1;
                if  assessConfig(counter, 3) == 0
                    assessConfig(counter, 1:2) = sysparam.offset;
                    assessConfig(counter, 3) = sysparam.tolerance;
                end
                
                summaryMatrix{counter} = [summaryMatrix{counter} metric.summaryMtx];
            end
        end
    end

%     summary = sum(summaryMatrix');


    % generate the export data
    if exist(exportPath, 'file')
        % if it exists already
        fid = fopen(exportPath, 'a');
    else
        fid = fopen(exportPath, 'a');

        % make headers
        header = ['offsetVal,offsetVal,errorTolerance,' ...
            'total,correct,falsePos,falseNeg,falsePos_OutOfSet,falseNeg_toofar,falseNeg_missing,TPPercentage,F1Seg'];

        fprintf(fid, [header '\n']);
    end

    for i = 1:length(summaryMatrix)
        %         fprintf(fid, [datestr(now) ',' num2str(observation{i}.subject) ',' num2str(observation{i}.session) ',' observation{i}.dataFolder ',' ]);
        %         fprintf(fid, [num2str(observationDataSource{i}.time{1}(end)) ',' num2str(timeKeep(1)) ',' num2str(timeKeep(2)) ',' num2str(timeKeep(4)) ',']);
        %         fprintf(fid, '%2.5f,%2.5f,%2.5f,', sysparam.comparison.offset, sysparam.comparison.tolerance);
        %         fprintf(fid, '%d,', summaryMtx);

        summary = sum(summaryMatrix{i}, 2);
        offset = assessConfig(i, 1:2);
        tolerance = assessConfig(i, 3);
        
        fScore_Seg(i) = fScoreCalculate(summary(2), 0, summary(3), summary(4), 1);
        
        fprintf(fid, '%2.5f,%2.5f,%2.5f,', offset, tolerance);
        fprintf(fid, '%d,', summary);
        fprintf(fid, '%f,%f,', summary(2)/summary(1), fScore_Seg(i));
        
        fprintf(fid, '\n');
    end
    
    maxFScore = max(fScore_Seg)

    fclose(fid);

end