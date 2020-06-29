function calcTrajTrig
    % load data and calculate joint angles via trig
    close all

    plotStuff = 0;
    fillThres = 2*pi;
    dt = 0.01;
    
    dataSource = 'healthy1';
    outputBase = 'C:\Documents\MATLABResults\DataPlots\';
    outputBasePath = [outputBase dataSource '\jointOutput\'];
    
    specStruct.dataset = dataSource;
    specStruct.patient = [18];
    specStruct.session = [];    

%     specStruct.exerciseAcceptPrefix = {'HAAO_STD', 'HEFO_STD', 'KEFO_SIT', 'KHEF_STD', 'KHEF_SUP', ...
%         'HFEO_STD', 'KEFO_SUP', 'KFEO_STD', 'KFEO_SUP', 'SQUA_STD'}; % set used for paper    
%     specStruct.exerciseAcceptPrefix = {'HEFO_STD', 'KEFO_SIT', 'KHEF_STD', 'KHEF_SUP', ...
%         'HFEO_STD', 'KEFO_SUP', 'KFEO_STD', 'KFEO_SUP'}; % hip first, sagittal only    
%     specStruct.exerciseAcceptPrefix = {'SQUA_STD'}; % ankle first, sagittal only

%     specStruct.exerciseAcceptPrefix = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP', 'SQUA_STD', 'STSO_SIT'};
%     specStruct.exerciseAcceptPrefix = {'HFEO_SUP', 'KEFO_SIT', 'KHEF_SUP'}; kinChain = 'hipToAnkle';
%     specStruct.exerciseAcceptPrefix = {'SQUA_STD', 'STSO_SIT'}; kinChain = 'ankleToHip';

	specStruct.exerciseAcceptPrefix = {'kefo_sit_slo'}; kinChain = 'hipToAnkle';

    
    basePath = 'C:\Documents\aslab\data';
    switch lower(dataSource)
        case 'healthy1'
            pathToRawData = [basePath '\Lowerbody_healthy1_2011-11'];
            
        case 'healthy2'
            pathToRawData = [basePath '\Lowerbody_healthy2_2013-07']; % path to the raw data
            
        case 'tri1'
            pathToRawData = [basePath  '\Lowerbody_TRI1_2012-10'];
            
        case 'stjoseph1'
            pathToRawData = [basePath  '\Lowerbody_StJoseph1_2013-02\'];
    end


%     options.headerMocap = []; 
%     options.headerImu = [];
    % options.headerEkf = []; 
    % options.headerSegManual = [];
    % options.headerSegCrop = [];
    options.headerSegAlg = [];
   
    fileStack = loadPatientFilepaths(pathToRawData, specStruct);
    
    % make sure the output filepath exist
    checkMkdir(outputBasePath);

    for ind_subjectData = 1:length(fileStack)
        % load each file separately so the memory usage isn't insane
        currFile = fileStack{ind_subjectData};
        
        savePath = fullfile(currFile.filePath, 'EKF', '2016_11_14_Planar');
        checkMkdir(savePath);
        saveTargetFile = fullfile(savePath, ['jointAngles.mat']);
        
        % load this specific EKF path
        %         options.headerEkf = fullfile(currFile.filePath, 'EKF', '2013_06_13', ['JointAngle_' currExerciseName '.data']);
        options.headerEkf = fullfile(currFile.filePath, 'EKF', '2015_03_23', 'ekf.header');
        
        %         options.headerSegManual = fullfile(currFile.filePath, 'Segmentation_manual_annotatedZVC',  'SegmentData_Manual_Manual.header');
        
        fprintf('%s (%u/%u): Reading %s\n', datestr(now), ind_subjectData, length(fileStack), currFile.filePath);
        subjectData = exerciseDataHandle(currFile.filePath, options);
        
        if isempty(subjectData)
            fprintf('  Skipping, no data loaded\n');
            continue
        end
        
         subjectData.filepathHeaderSegManual = fullfile(subjectData.dirPathExercise, 'Segmentation_manual_annotatedZVC',  'SegmentData_Manual_Manual.header');
        
        if ~exist(subjectData.filepathHeaderSegManual, 'file') && strcmpi(dataSource, 'healthy1')
            % no manual segments exists at all for this subject
            subjectData.filepathHeaderSegManual = fullfile(subjectData.dirPathExercise, 'Segmentation_manual',  'SegmentData_Manual_Manual.header');
        end
            
        if ~exist(subjectData.filepathHeaderSegManual, 'file')
            % no manual segments exists at all for this subject
            subjectData.filepathHeaderSegManual = [];
        end
        
        if isempty(subjectData.filepathHeaderSegManual) || ...
                ~exist(subjectData.filepathHeaderSegManual, 'file') || ...
                ~exist(subjectData.filepathHeaderEkf, 'file')
            fprintf('  Skipping, missing EKF or manseg\n');
            continue
            %             subjectData.filepathHeaderSegManual = manSegPath1;
        end
        
        subjectData.load;
        
        subjStr = num2str(subjectData.subjectName);
        if length(subjStr) == 1
            subjStr = ['0' subjStr];
        end
        
        sessStr = num2str(subjectData.session);
        if length(sessStr) == 1
            sessStr = ['0' sessStr];
        end
        
        %         subjectBlurb = [dataSource '_' ...
        %             'Subj' subjStr '_' ...
        %             'Sess' sessStr '_' ...
        %             subjectData.exerciseName];
        
        subjectBlurb = [upper(dataSource) '_' ...
            subjectData.exerciseType '_' ...
            'Subj' subjStr '_' ...
            'Sess' sessStr '_' ...
            subjectData.exerciseName];
        
        if isempty(subjectData.dataEkf) %|| isempty(subjectData.dataSegManual) %|| isempty(subjectData.dataSegCrop)
            
        else
            % DO SOMETHING
            
            marker = {};
            markerNew = {};
            joint = [];
            
            startTimeVal = mean([subjectData.dataSegCrop.timeStart]);
            endTimeVal = mean([subjectData.dataSegCrop.timeEnd]);
            
            [marker, markerName, trcTime, trcData, jointTime, jointData] = cleanMarkers(subjectData, startTimeVal, endTimeVal, dt, kinChain);
            
            % check for marker problems...if exist, remove the edge ones
            markersToCheck = abs([marker{:}]);
            markerCombined = sum(markersToCheck, 2);
            if ~isempty(find(markerCombined > 1e4)) || ~isempty(find(markerCombined  == 0))
                fprintf('  Skipping first and last segment\n');
                startTimeVal = mean([subjectData.dataSegManual.timeEnd(1) subjectData.dataSegManual.timeStart(2)]);
                endTimeVal = mean([subjectData.dataSegManual.timeEnd(end-1) subjectData.dataSegManual.timeStart(end)]);
                [marker, markerName, trcTime, trcData, jointTime, jointData] = cleanMarkers(subjectData, startTimeVal, endTimeVal, dt, kinChain);    
            end
            
% % %             % check for overage and crop along the segment edges
% % %             markerCombined = max([marker{:}], [], 2);
% % %             segmentsToKeep = 1:length(subjectData.dataSegManual.timeStart);
% % %             allExceed = find(markerCombined > 1e3);
% % %             for ind_exceeded = 1:length(allExceed)
% % %                 % find the corresponding segment block
% % %                 for ind_seg = 1:length(subjectData.dataSegManual.timeStart)
% % %                     if t
% % %                 end
% % %             end
            
            
            
            rotMtxLocal = rotx(1*pi/2);
%             rotMtxLocal = eye(3);
            switch subjectData.subjectName
                case {1, 19}
                    rotMtxLocal = rotMtxLocal*rotz(pi/2); % rotz -pi/2
            end
            
            for ind_la = 1:length(marker)
                marker{ind_la} = (rotMtxLocal*marker{ind_la}')';
            end
            
            T = viewmtx(0,0);
            T = T(1:3, 1:3);
            
            if plotStuff
                h1 = figure('Position', [   306   264   560   420]);
                h2 = figure('Position', [   803   264   560   420]);
            end
            
            for t = 1:1:length(marker{1}) % [7409 8772] %
                for i = 1:length(marker)
                    currMarker = marker{i}(t, :)';
                    currMapping(i, :) = T*currMarker(1:3);
                    markerNew{i}(t, :) = currMapping(i, [1 3]);
                end
                
                markerArray = [];
                for i = 1:length(marker)
                    markerArray = [markerArray; markerNew{i}(t, :)];
                end
                
                if plotStuff
                    figure(h1);
                    clf;
                    
                    hold on;
                    grid on;
                    
                    if 0
                        scatter3(currMapping(:, 1), currMapping(:, 2), currMapping(:, 3));
                        plot3(currMapping(:, 1), currMapping(:, 2), currMapping(:, 3));
                        scatter3(currMapping(1, 1), currMapping(1, 2), currMapping(1, 3), 'g');
                        
                        xlabel('x');
                        ylabel('y');
                        zlabel('z');
                        
                        view([160 15]);
                        
                        xlim([-1 1]);
                        ylim([-1 1]);
                        zlim([-1 1]);
                    else
                        scatter(markerArray(:, 1), markerArray(:, 2));
                        scatter(markerArray(1, 1), markerArray(1, 2), 'g');
                        plot(markerArray(:, 1), markerArray(:, 2));
                        
                        xlabel('x');
                        ylabel('y');
                        
                        xlim([-0.5 1.0]);
                        ylim([-1.0 0.5]);
                    end
                     
                    title([num2str(t) '/' num2str(length(marker{1}))]);
                    
                    figure(h2);
                    clf;
                    
                    plot(jointTime, jointData); hold on
                    [closeVal, closeInd] = findClosestValue(trcTime(t), jointTime);
                    plot([closeVal closeVal], ylim, 'k');
                    xlim([closeVal-1 closeVal+1]);
                    
                    pause(0.05);
                end
            end
            
            joint(:, 1) = zeros(size(markerNew{1}(:, 1)));
            
            % markerOffset = zeros(size(trcData.data.RHEEL));
            % marker{1} = trcData.data.RHEEL - markerOffset; % floor, marker 0
            % marker{2} = trcData.data.RANK - markerOffset;
            % marker{3} = trcData.data.RKNEEO - markerOffset;
            % marker{4} = trcData.data.RHIP - markerOffset;
            % marker{5} = trcData.data.STERNLOW - markerOffset;
            % marker{6} = trcData.data.STERNLOW - markerOffset;
            % marker{7} = trcData.data.STERNUP - markerOffset;
            % marker{8} = trcData.data.RELBL - markerOffset;
            
            for i = 2:length(marker)
                if i == 1
                    marker_first = zeros(size(markerNew{i}));
                    marker_second = markerNew{i};
                else
                    marker_first = markerNew{i-1};
                    marker_second = markerNew{i};
                end
                
                %     joint(:, i) = calcAngle_crossprod(marker_first, marker_second);
                jointTemp = calcAngle_tan(marker_first, marker_second)';
                
                switch kinChain
                    case 'ankleToHip'
                        angOffset = 0;
                        
                    case 'hipToAnkle'
                        angOffset = 0;
                end
                
                switch i
                    case 2
                        
                    case 3
                        joint(:, i-1) = jointTemp;
                    case 6
                        switch subjectBlurb
                            case 'HEALTHY1_KEFO_SIT_SLO_Subj18_Sess01_KEFO_SIT_SLO1'
                                joint(:, i-1) = jointTemp + pi - sum(joint(:, 1:i-2), 2);
                            otherwise
                                joint(:, i-1) = jointTemp - angOffset - sum(joint(:, 1:i-2), 2);
                        end
                       
                    otherwise
                        joint(:, i-1) = jointTemp - angOffset - sum(joint(:, 1:i-2), 2);
                end
% % %                 if i == 2
% % %                     joint(:, i-1) = jointTemp;
% % %                 else
% % %                     joint(:, i-1) = jointTemp - sum(joint(:, 1:i-2), 2);
% % %                 end
            end
            
            dataAngles = unwrap(joint);
%             dataAngles = (joint);
            fillInRegion = [];
            
            for i = 2:size(dataAngles, 2)
                for j = 2:size(dataAngles, 1)
                    if abs(dataAngles(j, i) - dataAngles(j-1, i)) > fillThres %% big jump occured
                        dataAngles(j, i) = dataAngles(j-1, i);
                        fillInRegion = [fillInRegion j];
                    end
                end
            end
            
            fillInRegion = unique(fillInRegion);
            
            markerArray_plot = [marker{:}];
            markerNewArray_plot = [markerNew{:}];
            
            saveTargetFilePng1 = fullfile(outputBasePath, ['1_' subjectBlurb '.png']);
            saveTargetFileFig1 = fullfile(outputBasePath, ['1_' subjectBlurb '.fig']);
            saveTargetFilePng2 = fullfile(outputBasePath, ['2_' subjectBlurb '.png']);
            saveTargetFileFig2 = fullfile(outputBasePath, ['2_' subjectBlurb '.fig']);
            saveTargetFilePng3 = fullfile(outputBasePath, ['3_' subjectBlurb '.png']);
            saveTargetFileFig3 = fullfile(outputBasePath, ['3_' subjectBlurb '.fig']);
            
            h1 = figure;
            % plot(unwrap(joint), '.'); hold on
            subplot(211);
            plot(markerArray_plot);
            title('marker array');
            ylim([-1.5 1.5]);
            subplot(212);
            plot(markerNewArray_plot);
            title('marker array new');
            ylim([-1.5 1.5]);
            saveas(h1, saveTargetFilePng1);
            saveas(h1, saveTargetFileFig1);
            
            edgeInd = 6 ;
            
            h2 = figure; 
            plot(dataAngles, 'b'); hold on;
            plot(dataAngles(:, :), '.');
            plot(fillInRegion,  edgeInd*ones(size(fillInRegion)), 'r.');
            plot(fillInRegion, -edgeInd*ones(size(fillInRegion)), 'r.');
            % plot(dataAngles);
            title(subjectBlurb);
            ylim([-edgeInd edgeInd]);
            saveas(h2, saveTargetFilePng2);
            saveas(h2, saveTargetFileFig2);
           
            h3 = figure;
            for ind_lala = 1:4
                subplot(2, 2, ind_lala);
                plot(marker{ind_lala+1});
                title(markerName{ind_lala+1});
                ylim([-1.7 1.7]);
            end
            saveas(h3, saveTargetFilePng3);
            saveas(h3, saveTargetFileFig3);
            
            close all
            
            kinematicData.time = trcTime;
            kinematicData.joint = dataAngles;
            kinematicData.marker = markerNew;
            save(saveTargetFile, 'kinematicData');
        end
    end
end

function rot = rotMtx(ang)
    rot = [cos(ang) -sin(ang); 
        sin(ang) cos(ang)];
end

function ang = calcAngle_tan(marker1, marker2)
    for i = 1:size(marker1, 1)
        marker1_curr = [marker1(i, :)];
        marker2_curr = [marker2(i, :)];
        
        x = marker2_curr(:, 1) - marker1_curr(:, 1);
        y = marker2_curr(:, 2) - marker1_curr(:, 2);
%         ang(:, i) = atan(y/x);
       ang(:, i) = atan2(y, x);
        
        if isnan(ang(:, i))
            ang(:, i) = 0;
        end
    end
end

function ang = calcAngle_crossprod(marker1, marker2)
    for i = 1:size(marker1, 1)
        marker1_curr = [marker1(i, :) 0];
        marker2_curr = [marker2(i, :) 0];
        
        crs = cross(marker1_curr, marker2_curr);
        mag = norm(marker1_curr) * norm(marker2_curr);
        
        if mag > 0
            arg = crs/mag;
            ang(:, i) = asin(arg(3));
        else
            ang(:, i) = 0;
        end
    end
    
    ang = ang';
end

function marker = combineMarkers(data, name, trcTimeOrig, trcTimeInterp)
    markerraw = [data.([name '_x']) data.([name '_y']) data.([name '_z']); ];
    marker = interp1(trcTimeOrig, markerraw, trcTimeInterp);
end

function marker = mergeMarkers(data1, data2)
    marker = zeros(size(data1));
    for i = 1:size(data1, 1)
        for j = 1:size(data1, 2)
            marker(i, j) = mean([data1(i, j) data2(i, j)]);
        end
    end
end

function [marker, markerName, trcTime, trcData, jointTime, jointData] = cleanMarkers(subjectData, startTimeVal, endTimeVal, dt, kinChain)
[~, trcStartInd] = findClosestValue(startTimeVal, subjectData.dataMocap.time);
            [~, trcEndInd] = findClosestValue(endTimeVal, subjectData.dataMocap.time);
%             trcInd = trcStartInd:trcEndInd;
            trcTime_load = subjectData.dataMocap.time;
            trcTime = trcTime_load(trcStartInd):dt:trcTime_load(trcEndInd);
            trcData = subjectData.dataMocap; 
            
            [~, ekfStartInd] = findClosestValue(startTimeVal, subjectData.dataEkf.time);
            [~, ekfEndInd] = findClosestValue(endTimeVal, subjectData.dataEkf.time);
            ekfInd = ekfStartInd:ekfEndInd;            
            jointTime = subjectData.dataEkf.time(ekfInd);
            jointData = subjectData.dataEkf.Q(ekfInd, :);
            
            if 0
                figure;
                plot([subjectData.dataMocap.time(1) subjectData.dataMocap.time(end)], [0.5 0.5], 'b'); hold on
                plot([subjectData.dataEkf.time(1) subjectData.dataEkf.time(end)], [-0.5 -0.5], 'r');
                plot([subjectData.dataSegCrop.timeStart(1) subjectData.dataSegCrop.timeStart(1)], ylim, 'k');
                plot([subjectData.dataSegCrop.timeEnd(1) subjectData.dataSegCrop.timeEnd(1)], ylim, 'k');
                ylim([-1 1]);
            end
            
            marker_toe = combineMarkers(trcData.data, 'R_Toe', trcTime_load, trcTime);
            
            marker_ankle_lat = combineMarkers(trcData.data, 'R_Ankle_Lateral', trcTime_load, trcTime);
%             marker_ankle_med = combineMarkers(trcData.data, 'R_Ankle_Medial', trcInd);
%             marker_ankle = mergeMarkers(marker_ankle_lat, marker_ankle_med);
            marker_ankle = marker_ankle_lat;
            
            marker_knee_lat = combineMarkers(trcData.data, 'R_Knee_Lateral', trcTime_load, trcTime);
%             marker_knee_med = combineMarkers(trcData.data, 'R_Knee_Medial', trcInd);
%             marker_knee = mergeMarkers(marker_knee_lat, marker_knee_med);
            marker_knee = marker_knee_lat;
            
            marker_hip = combineMarkers(trcData.data, 'R_ASIS', trcTime_load, trcTime);
            marker_shoulder = combineMarkers(trcData.data, 'R_Shoulder', trcTime_load, trcTime);
            
            switch kinChain
                case 'ankleToHip'
                    markerOffset = (marker_toe);
                    marker{1} = (marker_toe - markerOffset)/1000;
                    marker{2} = (marker_ankle - markerOffset)/1000;
                    marker{3} = (marker_knee - markerOffset)/1000;
                    marker{4} = (marker_hip - markerOffset)/1000;
                    marker{5} = (marker_shoulder - markerOffset)/1000;
                    marker{6} = (marker_shoulder - markerOffset)/1000;
                    marker{7} = (marker_shoulder - markerOffset)/1000;
                    marker{8} = (marker_shoulder - markerOffset)/1000;
                    markerName{1} = 'marker_toe';
                    markerName{2} = 'marker_ankle';
                    markerName{3} = 'marker_knee';
                    markerName{4} = 'marker_hip';
                    markerName{5} = 'marker_shoulder';
                    markerName{6} = 'marker_shoulder';
                    markerName{7} = 'marker_shoulder';
                    markerName{8} = 'marker_shoulder';
                    
                case 'hipToAnkle'
                    markerOffset = (marker_shoulder);
                    marker{1} = (marker_shoulder - markerOffset)/1000;
                    marker{2} = (marker_hip - markerOffset)/1000;
                    marker{3} = (marker_knee - markerOffset)/1000;
                    marker{4} = (marker_ankle - markerOffset)/1000;
                    marker{5} = (marker_toe - markerOffset)/1000;
                    marker{6} = (marker_toe - markerOffset)/1000;
                    marker{7} = (marker_toe - markerOffset)/1000;
                    marker{8} = (marker_toe - markerOffset)/1000;
                    markerName{1} = 'marker_shoulder';
                    markerName{2} = 'marker_hip';
                    markerName{3} = 'marker_knee';
                    markerName{4} = 'marker_ankle';
                    markerName{5} = 'marker_toe';
                    markerName{6} = 'marker_toe';
                    markerName{7} = 'marker_toe';
                    markerName{8} = 'marker_toe';
            end
end