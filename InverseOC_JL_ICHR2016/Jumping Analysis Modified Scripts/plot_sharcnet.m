function al
clearvars
clc

    addpath(genpath(fullfile('..','Symoro')));
    addpath(genpath(fullfile('..','Model')));
    addpath(genpath(fullfile('..','Common')));
    addpath(genpath(fullfile('..','..','Toolboxes')));
    addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\MarkerSwappingEKF\2018_07_04')); % what folders are correct???
    addpath(genpath('C:\Users\kgwester\Documents\ResearchWork\MarkerSwappingEKF\ekf_2018_05_29')); % what folders are correct???

baseSourceFolder = 'C:\Users\kgwester\Documents\ResearchWork\ioc_project\';
commentsArray = {''};
SVDArray = {'10000'};
thresholdMultiplierArray = [1];

nowStrMaster = datestr(now, 'yyyymmddHHMMSS');

for ind_multiplierarray = 1:length(thresholdMultiplierArray)
for ind_commentsarray = 1:length(commentsArray)
    for ind_svdarray = 1:length(SVDArray)
        
        nowStr = datestr(now, 'yyyymmddHHMMSS');
        comments = commentsArray{ind_commentsarray};
        svdEntry = SVDArray{ind_svdarray};
        thresholdMultiplier = thresholdMultiplierArray(ind_multiplierarray);
        
basePath_array =           {};
basePathPlot_array =       {};
basePathMasterPlot_array = {};        
        
basePath_array{1} =           [];
basePathPlot_array{1} =       [];
basePathMasterPlot_array{1} = [];

basePath_array{end+1} =           [baseSourceFolder 'jump_results\jump2d\'];
basePathPlot_array{end+1} =       [baseSourceFolder 'jump_results\' 'p_' nowStr '_' comments '_' svdEntry '\'];
basePathMasterPlot_array{end+1} = basePathPlot_array{end};

for i = 1:length(basePath_array)
    if isempty(basePath_array{i})
        continue
    end
    
basePath = basePath_array{i};
basePathPlot = basePathPlot_array{i};
basePathMasterPlot = basePathMasterPlot_array{i};

masterSource = [basePathPlot 'log_overall.csv'];
masterFigureSave = [basePathPlot '\fig\'];

checkMkdir(basePathMasterPlot);
checkMkdir(masterFigureSave);

dir_layer1 = dir(basePath);

counter = 0;
totalIndWeighting = [];
accLess = {};
accMore = {};
maxIndTrack = [];
maxAccIndiv = [];
instanceNameCounter = {};
% rmse_reportCounter = [];
windowedRmse_mean = [];
windowedRmse_std = [];
windowedResnorm_mean = [];
windowedResnorm_std = [];
                    
currDataset = strsplit(basePath_array{i}, '\'); 
masterSaveLog = fullfile(basePathMasterPlot, ['outputFile_plot_' currDataset{end-2} '_' nowStr '_' comments '_' svdEntry '.csv']);
masterSaveLogSeg = fullfile(basePathMasterPlot, ['outputFile_segment_' currDataset{end-2} '_' nowStr '_' comments '_' svdEntry '.csv']);

masterSaveOverall = fullfile(baseSourceFolder, ['overallFile_' nowStrMaster '.csv']);
masterSaveOverallSeg = fullfile(baseSourceFolder, ['overallFile_segment_' nowStrMaster '.csv']);

for ind_layer1 = 1:length(dir_layer1)
    currDir_layer1 = dir_layer1(ind_layer1);
    currDir_layer1_fullName = [basePath currDir_layer1.name '\'];
    
    splitTheStr = strsplit(currDir_layer1.name, '_');
    
    if ~currDir_layer1.isdir
        continue
    elseif length(currDir_layer1.name) < 3
        continue
%     elseif ~strcmpi(splitTheStr{2}, '3r')
%         continue
    end
    
    if length(splitTheStr) > 8 && ~strcmpi(splitTheStr{9}, svdEntry)
        continue
    end
    
    dir_layer2 = dir(currDir_layer1_fullName);
    for ind_layer2 = 1:length(dir_layer2)
        currDir_layer2 = dir_layer2(ind_layer2);
        currDir_layer2_fullName = [currDir_layer1_fullName currDir_layer2.name '\'];
        
        if strcmpi(currDir_layer2.name, 'log_overall.csv')
            % copy it to the master source
            fId1 = fopen([currDir_layer2_fullName(1:end-1)]);
            fId2 = fopen(masterSource, 'a+');
            
            tline = fgetl(fId1);
            while ischar(tline)
                fprintf(fId2, '%s\n', tline);
                tline = fgetl(fId1);
            end
            
            fclose(fId1);
            fclose(fId2);
            
            continue
        elseif ~currDir_layer2.isdir
            continue
        elseif length(currDir_layer2.name) < 3
            continue
        end
        
        dir_layer3 = dir([currDir_layer2_fullName]);
        
        for ind_layer3 = 1:length(dir_layer3)
            currDir_layer3 = dir_layer3(ind_layer3);
            currDir_layer3_fullName = [currDir_layer2_fullName currDir_layer3.name '\'];
                
            if ~currDir_layer3.isdir
                continue
            elseif length(currDir_layer3.name) < 3
                continue
            end
            
            fprintf('Loading %s\n', currDir_layer3_fullName);
            
            dir_layer4 = dir([currDir_layer3_fullName '*end.mat']);
            matFilesToLoad = {};
            indStart = [];
            
            if isempty(dir_layer4)
                continue
            end
            
% % %             try
                for ind_layer4 = 1:length(dir_layer4)
                    currDir_layer4 = dir_layer4(ind_layer4);
                    matFilesToLoad{ind_layer4} = [currDir_layer3_fullName currDir_layer4.name];
                    
                    subtextString = strsplit(currDir_layer4.name, '_');
                    indStart(ind_layer4) = str2num(subtextString{4});
                    indEnd(ind_layer4) = str2num(subtextString{5});
                end
                
                [Y, I] = sort(indStart);
                indStart = indStart(I);
                indEnd = indEnd(I);
                indTotal = [indStart' indEnd'];
                matFilesToLoadUse = matFilesToLoad(I);
                
                instanceName = currDir_layer3.name;
%                 currDir_layer4_plotFullName = [basePathPlot currDir_layer1.name '\' currDir_layer2.name '\' currDir_layer3.name];
                currDir_layer4_plotFullName = [basePathPlot currDir_layer1.name '\' currDir_layer3.name];
           
                if ~isempty(matFilesToLoadUse)                
                    masterSave1 = fullfile(basePathMasterPlot, [currDir_layer1.name '_fig_' nowStr]);
                    loadFileAndPlot_overlapMode(matFilesToLoadUse, instanceName, indTotal, currDir_layer4_plotFullName, masterSaveOverall, masterSave1);
                end
                
%                 if ~isempty(matFilesToLoadUse)
%                     counter = counter + 1;
%                     instanceNameCounter{counter} = instanceName(end-12:end-1);
%                     [totalIndWeighting(counter), accLess{counter}, accMore{counter}, cost_function_names, constThres, maxIndTrack(counter), rmse_reportCounter] ...
%                         = loadFileAndSegment(matFilesToLoadUse, instanceName, indTotal, currDir_layer4_plotFullName, masterSaveLogSeg, masterSaveOverallSeg, comments, thresholdMultiplier, masterFigureSave);
%                     
%                     maxAccStructIndiv(counter) = rmse_reportCounter.accMax(rmse_reportCounter.dofMax);
%                     windowedRmse_mean(counter) = rmse_reportCounter.windowed_rmse_belowThreshold_mean;
%                     windowedRmse_std(counter) = rmse_reportCounter.windowed_rmse_belowThreshold_std;
%                     windowedResnorm_mean(counter) = rmse_reportCounter.resnorm_mean;
%                     windowedResnorm_std(counter) = rmse_reportCounter.resnorm_std;
%                 end
                
                fclose all;
        end
    end
end
end
end
end
end
end
