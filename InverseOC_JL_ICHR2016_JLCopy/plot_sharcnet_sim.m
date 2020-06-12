% look into all folders and plot the results. designed to plot the sharcnet
% results
clearvars
clc

baseSourceFolder = 'D:\results\IOC\IOC10_modelRL2D_qVB\';
nowStr = datestr(now, 'yyyymmddHHMMSS');
    
basePath_array =           {};
basePathPlot_array =       {};
basePathMasterPlot_array = {};        
        
basePath_array{1} =           [];
basePathPlot_array{1} =       [];
basePathMasterPlot_array{1} = [];

basePath_array{end+1} =           [baseSourceFolder 'sim\CF3\'];
basePathPlot_array{end+1} =       [baseSourceFolder 'sim\CF3' 'p_' nowStr '\'];
basePathMasterPlot_array{end+1} = basePathPlot_array{end};


for i = 1:length(basePath_array)
    if isempty(basePath_array{i})
        continue
    end
    
basePath = basePath_array{i};
basePathPlot = basePathPlot_array{i};
basePathMasterPlot = basePathMasterPlot_array{i};

checkMkdir(basePathMasterPlot);

masterSource = [basePathPlot 'log_overall.csv'];

dir_layer1 = dir(basePath);

counter = 0;
totalIndWeighting = {};
accLess = {};
accMore = {};

currDataset = strsplit(basePath_array{i}, '\'); 
masterSaveLog = fullfile(basePathMasterPlot, ['outputFile_plot_' currDataset{end-2} '_' nowStr '.csv']);
masterSaveLogSeg = fullfile(basePathMasterPlot, ['outputFile_segment_' currDataset{end-2} '_' nowStr '.csv'])

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
                    indStart(ind_layer4) = str2num(subtextString{6});
                    indEnd(ind_layer4) = str2num(subtextString{7});
                end
                
                [Y, I] = sort(indStart);
                indStart = indStart(I);
                indEnd = indEnd(I);
                indTotal = [indStart' indEnd'];
                matFilesToLoadUse = matFilesToLoad(I);
                
                instanceName = currDir_layer3.name;
%                 currDir_layer4_plotFullName = [basePathPlot currDir_layer1.name '\' currDir_layer2.name '\' currDir_layer3.name];
                currDir_layer4_plotFullName = [basePathPlot currDir_layer1.name '\' currDir_layer3.name];
             
%                 if strcmpi(currDir_layer3.name(1:5), 'Subj1') && strcmpi(currDir_layer1.name(1:7), 'win2011')
%                     continue
%                 end

%                 if ~strcmpi(currDir_layer3.name(1:5), 'Subj6')
%                     continue
%                 end

                if ~isempty(matFilesToLoadUse)                
                    masterSave1 = fullfile(basePathMasterPlot, [currDir_layer1.name '_fig_' nowStr]);
                    loadFileAndPlot_overlapMode(matFilesToLoadUse, instanceName, indTotal, currDir_layer4_plotFullName, masterSaveLog, masterSave1);
                end
                
%                 if ~isempty(matFilesToLoadUse)
%                     counter = counter + 1;
%                     
%                      [totalIndWeighting{counter}, accLess{counter}, accMore{counter}, cost_function_names] ...
%                         = loadFileAndSegment(matFilesToLoadUse, instanceName, indTotal, currDir_layer4_plotFullName, masterSaveOverallSeg, comments);
%                 end
                
                fclose all;
                
                %             tallyfileName = [currDir_layer2_fullName instance Name '_tally.csv'];
                %             targetFileName = [currDir_layer3_plotFullName '\' instanceName '_tally.csv'];
                %
                %             if ~strcmpi(tallyfileName, targetFileName)
                %                 copyfile(tallyfileName, targetFileName);
                %             end
% % %             catch err
% % %                 close all;
% % %                 err.message
% % %             end
            
            % copy the tally file over
        end
    end
end

end
