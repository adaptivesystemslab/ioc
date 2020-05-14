function calcTrajTrig
% load data and calculate joint angles via trig

plotStuff = 1;  
fillThres = 0.7;

close all
marker = {};
markerNew = {};
joint = [];

masterPath = ['C:\Documents\aslab\data\APARS_Exploratory\2011-02-07 LegRaise_MA'];
loadPath = [masterPath '\Jon_Rand1.trc'];
trcData_load = readTrc(loadPath);

% jointAngleFile = searchForFileByExt(fullfile(masterPath, 'JointAngles', 'IK_2016-01_MK'), 'sub*.mat');
% loadPath = fullfile(fullfile(masterPath, 'JointAngles', 'IK_2016-01_MK'), jointAngleFile);
% jointData_load = load(loadPath);
% jointData_load = jointData_load.dataAngles;
% dofMapping = [38 27 34 33 38 36 23];
% jointData_load = [jointData_load zeros(size(jointData_load, 1), 1)];
% jointData = jointData_load(:, dofMapping);

savePath = [masterPath '\JointAngles\IK_2016-10_JL\'];
saveTargetFile = [savePath 'jointAngles.mat'];
checkMkdir(savePath);

trcTime = trcData_load.data.Time;
trcData = trcData_load.data;

markerOffset = (trcData.Hip);
marker{1} = (trcData.Hip - markerOffset)/1000;
marker{2} = (trcData.L_Knee - markerOffset)/1000;
marker{3} = (trcData.L_Ankle - markerOffset)/1000;

rotMtxLocal = rotx(1*pi/2);
for ind_la = 1:length(marker)
    marker{ind_la} = (rotMtxLocal*marker{ind_la}')';
end

jointTime = trcData.Time;  
jointData = [marker{:}];

% chain for 2015 vs mapping in 2011
% 1. Foot      N/A (38)
% 2. Ankle     R FOOT Z (27)
% 3. Knee      R KNEE X (34)
% 4. Hip       R HIP Z (33)
% 5. Pelvis    N/A (38)
% 6. Torso     U BODY Y (36)
% 7. Shoulder  R SHOULDER Z (23)

% param.L1=norm(markers(23).translation(2)); % RANK (23) toe to ankle?
% param.L2=norm(markers(23).translation(1)); %
% param.L3=norm(markers(23).translation - markers(30).translation); % RANK - RKNEE
% param.L4=norm(markers(30).translation - markers(5).translation); % RKNEE - RHIP
% param.L5=0;  % RHIP - ASIS (but set to 0)
% param.L6=norm(markers(5).translation - markers(32).translation);  % ASIS to STN (but set to RHIP - STERNLOW)
% param.L7=norm(markers(32).translation - markers(31).translation); % STR to C7/CLA (but set to STERNLOW to STERNUP)
% param.L8=norm(markers(31).translation - mean([markers(21).translation; markers(22).translation])); % STERNUP - RELB

% markerOffset = zeros(size(trcData.data.RHEEL));
markerOffset = (trcData.Hip);
marker{1} = (trcData.Hip - markerOffset)/1000;
marker{2} = (trcData.L_Knee - markerOffset)/1000;
marker{3} = (trcData.L_Ankle - markerOffset)/1000;


rotMtxLocal = rotx(1*pi/2);
for ind_la = 1:length(marker)
    marker{ind_la} = (rotMtxLocal*marker{ind_la}')';
end

T = viewmtx(0,0,0);
            T = T(1:3, 1:3);
            
            if plotStuff
                h1 = figure('Position', [   306   264   560   420]);
                h2 = figure('Position', [   803   264   560   420]);
            end
            
            for t = 1:1:length(marker{1}) % [7409 8772] %
                for i = 1:length(marker)
                    currMarker = marker{i}(t, :)';
                    currMapping(i, :) = T*currMarker(1:3);
                    markerNew{i}(t, :) = currMapping(i, 1:2);
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
                    
                    if 1
                        scatter3(currMapping(:, 1), currMapping(:, 2), currMapping(:, 3));
                        plot3(currMapping(:, 1), currMapping(:, 2), currMapping(:, 3));
                        scatter3(currMapping(1, 1), currMapping(1, 2), currMapping(1, 3), 'g');
                        
                        xlabel('x');
                        ylabel('y');
                        zlabel('z');
                        
                        view([100 10]);
                        
                        xlim([-1 1]);
                        ylim([-1 1]);
                        zlim([-1 1]);
                    else
                        scatter(markerArray(:, 1), markerArray(:, 2));
                        scatter(markerArray(1, 1), markerArray(1, 2), 'g');
                        plot(markerArray(:, 1), markerArray(:, 2));
                        
                        xlabel('x');
                        ylabel('y');
                        
                        xlim([-0.5 0.5]);
                        ylim([-0.2 1.5]);
                    end
                    
                    title([num2str(t) '/' num2str(length(marker{1}))]);
                    
                    figure(h2);
                    clf;
                    
                    plot(jointTime, jointData); hold on
                    [closeVal] = findClosestValue(trcTime(t), jointTime);
                    plot([closeVal closeVal], ylim, 'k');
%                     xlim([closeVal-1 closeVal+1]);
                    
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
                if i == 2
                    joint(:, i-1) = jointTemp;
                else
                    joint(:, i-1) = jointTemp - sum(joint(:, 1:i-2), 2);
                end
                
%                 for j = 1:1:5000
%                     lenN = norm(marker_second(j, :) - marker_first(j, :));
%                     len = [lenN 0]';
%                     rot_curr = rotMtx(sum(joint(j, :)));
%                     ef_prev = ef{i-1}(j, :)';
%                     ef_curr(j, :) = ef_prev + rot_curr*len;
%                     ef{i} = ef_curr;
%                 end
            end
            
            % figure;
            % for t = 1:10:5000
            % % scatter(ef{1}(1, t), ef{1}(2, t), 'o');
            % for i = 3:length(ef)
            %     hold on
            %     scatter(ef{i}(1, t), ef{i}(2, t), 'o');
            %     plot([ef{i-1}(1, t) ef{i}(1, t)], [ef{i-1}(2, t) ef{i}(2, t)], '-');
            % end
            % end
            
            dataAngles = unwrap(joint);
%             dataAngles = (joint);
            fillInRegion = [];
            
            for i = 2:4 %1:size(dataAngles, 2)
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

h = figure;
% plot(unwrap(joint), '.'); hold on
plot(dataAngles, 'b'); hold on
plot(dataAngles(:, 2:4), '.');
plot(fillInRegion, -3*ones(size(fillInRegion)), 'r.');
ylim([-3 3]);

saveas(h, ['C:\Documents\MATLABResults\DataPlots\squats_tuat_2011\Subj' subjN '.png']);
save(saveTargetFile, 'dataAngles');
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