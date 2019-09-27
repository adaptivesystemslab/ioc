% Reads in previously-saved trc and ekf marker position data, replays
% visualization using ekf_state and lg_ekf_hatInvState (joint angles)

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
addpath('..\..\');
addpath(EKFCodePath);


partNum = '17';
i_targ = '85';
i_set = '1';
i_jump = '1';


makeMovie = 1;
movieName = ['P' partNum '_' i_targ '_' i_set '_' i_jump '_markers_only.avi'];

dataFolder = ['P' partNum '_filtered/'];
data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));



mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 
loadFilePath = [newdir '\results\RESULTS_EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '_const_model.mat'];

% if(exist(loadFilePath,'file')~=2)
%     disp(['EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '.mat NOT FOUND']);
% else
    load(loadFilePath);


% Read in .trc files and make models

dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P' partNum '_template'];
% dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P']; %because some files names "P03-template"
% file_ind = strfind(data_files,dataName);
% file_ind = find(not(cellfun('isempty', file_ind)));
% dataName = data_files{file_ind};

dataNameSL = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-start_land'];

%% Read in marker data, make models, visualize
trc = readTrc(['../data/' dataFolder dataName '.trc']);
markerNames = fieldnames(trc.data);
markerNames = markerNames(3:end); % first 2 names are Frame # and Time
% Rotate markers so subject jumps in positive X direction
for m = 1:length(markerNames)
    trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
end
trc = fillInTRCFrames(trc);


trcSL = readTrc(['../data/' dataFolder dataNameSL '.trc']);
markerNamesSL = fieldnames(trcSL.data);
markerNamesSL = markerNamesSL(3:end); % first 2 names are Frame # and Time
% Rotate markers so subject jumps in positive X direction
for m = 1:length(markerNamesSL)
    trcSL.data.(markerNamesSL{m}) = (rotz(-pi/2)*trcSL.data.(markerNamesSL{m})')';
end
trcSL = fillInTRCFrames(trcSL);



% Make Visualizer
vis = rlVisualizer('vis',640,960);
vis.update();


if(makeMovie)
    vidObj = VideoWriter(movieName);
    vidObj.Quality = 100;
    vidObj.FrameRate = 20; % 200Hz data, plots every 10 frames
    open(vidObj);
end



for i= 1:size(trc.data.Frame,1)

    if(mod(i,5)==0)
%             disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);

        z_eul = mes_eul(i,:)';

        % Display Markers
        for m = 1:numel(z_eul)/3
            m_pos = z_eul(m*3-2:m*3);
            vis.addMarker(num2str(m),m_pos); %,[1 1 0 1]);
        end
        vis.update();
        
        for m = 1:numel(markerNamesSL)
            m_pos = (trcSL.data.(markerNamesSL{m})(i,:))/1000;
            vis.addMarker([ 'G' num2str(m)],m_pos,[1 0 0 1]);
        end

        vis.update();

        if(makeMovie)
            pause(0.001);
            [img,imgbool] = vis.getScreenshot();
%                 img = flipdim(img,1);
            imshow(img);
            writeVideo(vidObj, getframe(gca)); 
        end
    end

end

clear vis;

if (makeMovie) 
    close(vidObj); 
end
