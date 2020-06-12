function [crossingStruct, crossingDof] = zvcCheck(timeData, veloToTest, crossingStruct, timeStep, threshold)
    % Check for zero velocity crossings
    
    W = length(veloToTest);
    
    if ~exist('crossingStruct', 'var') || isempty(crossingStruct)
        dof = size(veloToTest, 1);
        
        for i = 1:dof            
            crossingStruct{i}.Index = []; % TODO check this later. Start each crossing logger with a ZVC
            crossingStruct{i}.Time = [];
            crossingStruct{i}.Mtx = [];
            crossingStruct{i}.prevCrossingInd = -1;
        end
    end
    
    if ~exist('timeStep', 'var')
        timeStep = 0; %size(veloToTest, 2);
    end
    
    if ~exist('threshold', 'var')
        rangeFull = range(veloToTest);
        threshold = 0.01*rangeFull;
    end
    
    % window step for ZVC sweeping
%     zvcWindowStep = W-1; % this would be the full size of the W
    zvcWindowStep = 1;
    zvcWindowRange = 3;
    
    distanceFromPreviousZVC = 10; 
    
    crossingDof = zeros(size(veloToTest, 1), 1);
    iInd = 1:size(veloToTest, 1);
    jInd = 1:zvcWindowStep:size(veloToTest, 2)-1;
    
    for i = iInd % each dof...
%         range = max(veloToTest(i, 2:end-1)) - min(veloToTest(i, 2:end-1));
        crossing = zeros(1, length(jInd));

        for j = jInd % 2 timestep window       
            if j+zvcWindowRange < W
                rangeExamine = abs(mean(veloToTest(i, j:j+zvcWindowRange)));
                if rangeExamine < threshold
                    % Low velocity entry (change this to 1 and comment out the
                    % bottom stuff in order to switch back to previous mode)
                    crossing(j) = 1; % TODO change to the entire window?
                end
            end
            
            % actual ZVCs is more important
            if veloToTest(i, j) < 0 && veloToTest(i, j+zvcWindowStep) > 0
                % ZVC: if a positive crossing was made
                crossing(j) = 1;
            elseif veloToTest(i, j) > 0 && veloToTest(i, j+zvcWindowStep) < 0
                % ZVC: if a negative crossing was made
                crossing(j) = 1; % was -1 (direction mattered) before %TODO make a decision
            end

            % insert it into the crossing matrix, if it passes the time treshold
            if crossing(j) ~= 0 && timeStep + j >= (crossingStruct{i}.prevCrossingInd + distanceFromPreviousZVC)
                % positive threshold observed
%                 crossingStruct.Counter(i) = crossingStruct.Counter(i) + 1;
    
                nextZvcInd = timeStep + j + floor(zvcWindowStep/2);
                
                if nextZvcInd == 0
                    nextZvcInd = 1;
                elseif nextZvcInd > length(timeData)
                    nextZvcInd = length(timeData);
                end

                currCrossingInd = [crossingStruct{i}.Index nextZvcInd];
                crossingStruct{i}.Index = sort(unique(currCrossingInd)); % feed into timeData to get crossingTime
                crossingStruct{i}.prevCrossingInd = timeStep + j; % update the 'last crossing' item
                
                crossingDof(i) = 1; % the dof that got a crossing
            end
        end

        currCrossingInd = crossingStruct{i}.Index;
        crossingStruct{i}.Time = timeData(currCrossingInd);
        crossingStruct{i}.Mtx = ones(size(currCrossingInd));
        crossingStruct{i}.prevCrossingChar = crossing;
        crossingStruct{i}.Counter = length(currCrossingInd);
    end 
end