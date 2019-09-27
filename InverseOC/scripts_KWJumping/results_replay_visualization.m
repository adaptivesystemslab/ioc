% Reads in previously-saved trc and ekf marker position data, replays
% visualization using ekf_state and lg_ekf_hatInvState (joint angles)

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
addpath('..\..\');
addpath(EKFCodePath);

partNum = '07';
i_targ = '70';
i_set = '1';
i_jump = '3';

plotEveryNFrames = 5;
plotMarkers = 0;
use_const_model = 1;

dataFolder = ['P' partNum '_filtered/'];
data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));

makeMovie = 0;
movieName = ['P' partNum '_' i_targ '_' i_set '_' i_jump '_kin_model.avi'];



mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 
if(use_const_model)
    loadFilePath = [newdir '\results\RESULTS_EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '_const_model.mat'];
else
    loadFilePath = [newdir '\results\RESULTS_EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '.mat'];
end

if(exist(loadFilePath,'file')~=2)
    disp(['EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '.mat NOT FOUND']);
else
    load(loadFilePath);
    % this loads in: 'mes_eul','mes_eul_predict',...
    %    'mes_lie','mes_lie_predict','ekf_state','lg_ekf_hatInvState'
    
    [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    % Read in .trc files and make models
    
    dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P' partNum '_template'];
%     dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P']; %because some files names "P03-template"
%     file_ind = strfind(data_files,dataName);
%     file_ind = find(not(cellfun('isempty', file_ind)));
%     dataName = data_files{file_ind};
    
%     dataNameSL = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
%                 '_' num2str(i_jump) '_clean-start_land'];

    %% Read in marker data, make models, visualize
    trc = readTrc(['../data/' dataFolder dataName '.trc']);
    markerNames = fieldnames(trc.data);
    markerNames = markerNames(3:end); % first 2 names are Frame # and Time
    % Rotate markers so subject jumps in positive X direction
    for m = 1:length(markerNames)
        trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
    end
    trc = fillInTRCFrames(trc);
    
    
    %Create Euler and Lie Group Models
    initPosWindow = [251, 450];
%     [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_setLinkLengths(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath,partNum);
    [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_setLinkLengths_shldrPrism(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath,partNum);

    mdl_eul.forwardPosition;
    


    % Make Visualizer
    vis = rlVisualizer('vis',800,700);
    vis.addModel(mdl_eul);
%     vis.addModel(mdl_lie);
    vis.update();
    
    %Initialize model, set to initPos (approximate start position)
    mdl_eul.position = initPos;
    mdl_eul.velocity(:) = 0;
    mdl_eul.acceleration(:) = 0;
    mdl_eul.forwardPosition;
    mdl_eul.forwardVelocity;


    vis.update;

    ekf = EKF_Q_DQ_DDQ(mdl_eul);
    
    
    if(makeMovie)
        vidObj = VideoWriter(movieName);
        vidObj.Quality = 100;
        vidObj.FrameRate = 20; % 200Hz data, plots every 5 frames, 50% speed
        open(vidObj);
    end
    
    
    for i= 1:size(mes_eul,1)

        % Set EKF model state
        ekf.makeMeasure(ekf_state(i,:));
        %NOTE: "makeMeasure" function call includes model position update

        if(mod(i,plotEveryNFrames)==0)
%             disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
            
            if(plotMarkers)
                z_eul = mes_eul(i,:)';
%                 z_lie = mes_lie(i,:)';
                
                % Display Markers
                for m=1:numel(z_eul)/3
                    m_pos = z_eul(m*3-2:m*3);
                    vis.addMarker(num2str(m),m_pos,[1 1 0 1]);
                end

                
            end
            
            vis.update();
            
            if(makeMovie)
                if(i==plotEveryNFrames) % pause to set visualizer frame 
                    pause(0.001);
                end
                pause(0.001);
                [img,imgbool] = vis.getScreenshot();
    %                 img = flipdim(img,1);
                imshow(img);
                writeVideo(vidObj, getframe(gca)); 
            else
                pause(0.015);
            end
        end
        
    end
    
    clear vis;
    
    if (makeMovie) 
        close(vidObj); 
    end
end