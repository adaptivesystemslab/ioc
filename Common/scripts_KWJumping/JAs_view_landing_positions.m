% load JA data from multiple participants, get all jump grades from all
% jumps into single array



partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
targetDist = {'55','70','85'};
% ballOfFootMarkers = [21, 22, 28, 29]; % FLLat, FLMed, FRLat, FRMed
ballOfFootMarkers = {'FLLat', 'FLMed', 'FRLat', 'FRMed'};

jumpGradeTypes = {'B','SB','P','P*','SF','F'};
load('jumpGradesAllParts.mat');

colorVec = [0.1,0.1,1; 0.8,0.5,1; 0,0.9,0; 0,0.9,0; 0.9,0.7,0; 0.9,0.1,0];
%        = [dark blue, light purple, green, green, yellow/orange, red]


figure(1); clf; hold on; grid on;
axis([-0.25, 0.25, -0.3, 0.2]); 
plot([0, 0],[ -0.3, 0.2], 'k--');
plot([-0.25, 0.25],[0, 0], 'k--');
rectangle('Position',[-0.2, -0.1, 0.4, 0.2]);



for partNum = 1:numel(partsToRun)
    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partsToRun{partNum} '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partsToRun{partNum} '.mat NOT FOUND']);
    else
        load(loadFilePath); % JA data struct, human data
        
        dataFolder = ['P' partsToRun{partNum} '_filtered/'];
        data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));
        
        for i_targ = 1:3
            for i_set = 1:2
                for i_jump = 1:6
                    currJump = 12*(i_targ-1) + 6*(i_set-1) + i_jump;
                    
                    dataName = ['P' partsToRun{partNum} '_target_' targetDist{i_targ} '_' num2str(i_set) ...
                                        '_' num2str(i_jump) '_clean-P']; %because some files names "P03-template"
                    file_ind = strfind(data_files,dataName);
                    file_ind = find(not(cellfun('isempty', file_ind)));
                    dataName = data_files{file_ind};


                    trc = readTrc(['../data/' dataFolder dataName]);
                    markerNames = fieldnames(trc.data);
                    markerNames = markerNames(3:end); % first 2 names are Frame # and Time
                    % Rotate markers so subject jumps in positive X direction
                    for m = 1:length(markerNames)
                        trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
                    end
                    
                    
                    
                    currGrade = jumpGradesAll_num(currJump, partNum);
                    landFrame = JA.LandFrame(6*(i_set-1) + i_jump, i_targ);
                    landPos = JA.locationLand(currJump);
                    
                    footPos = zeros(4,1);
                    for m = 1:numel(ballOfFootMarkers)
                        footPos(m,1:2) = trc.data.(ballOfFootMarkers{m})(landFrame+20,1:2);
                    end
                    
                    footPosFwd = (mean(footPos(:,1))/1000) - landPos;
                    footPosSide = mean(footPos(:,2))/1000;
                    
                    if(currGrade==3 || currGrade==4)
                        plot(footPosSide, footPosFwd, '.', 'Color', colorVec(currGrade,:), 'MarkerSize', 5);
                    else
                        plot(footPosSide, footPosFwd, '.', 'Color', colorVec(currGrade,:), 'MarkerSize', 10);
                    end
                    
                end
            end
        end
        
        
        
    end
end

