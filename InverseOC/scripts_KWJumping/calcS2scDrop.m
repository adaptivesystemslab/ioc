%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
clear vis; close all;

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);

visualize = 0;
plotJointAngles = 0;
pauseActive = 0;
saveMarkerData = 1;
saveJointAngles = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partsToRun = {'02','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'04'};

dataFolder = 'Shoulder_Center_Recordings/';

data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));

S2sc_L = zeros(22,1);
S2sc_R = zeros(22,1);

% Large Nested FOR Loops for each jump recording
for i = 1:numel(partsToRun)
    partNum = partsToRun{i};
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = -S2sc*10; % convert to mm, since .TRC data is in mm
    
    dataName = ['P' partNum '_shoulder_center_clean.trc'];
    %             file_ind = strfind(data_files,dataName);
    %             file_ind = find(not(cellfun('isempty', file_ind)));
    %             dataName = data_files{file_ind};

    if(exist(['../data/' dataFolder dataName],'file')~=2)
        disp([dataName '.TRC does not exist'])

    else
        %% Read in marker data, make models, visualize
    %                 trc = readTrc(['../data/' dataFolder dataName '.trc']);
        trc = readTrc(['../data/' dataFolder dataName]);
        markerNames = fieldnames(trc.data);
        markerNames = markerNames(3:end); % first 2 names are Frame # and Time
        % Rotate markers so subject jumps in positive X direction
        for m = 1:length(markerNames)
            trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
        end

        
        window = [size(trc.data.Neck,1)/2,size(trc.data.Neck,1)];
        [S2sc_L(str2num(partNum)), S2sc_R(str2num(partNum))] = getS2scDrop(trc,0,window);
        
    end
end




