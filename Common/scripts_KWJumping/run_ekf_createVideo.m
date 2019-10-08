%Full Body LG EKF and EUL EKF for the boxing CMU data
%Load up the generated model
clear vis; close all;

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';

addpath(genpath(pwd));
addpath(EKFCodePath);

visualize = 1;
makeMovie = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partNum = '17';
[age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
S2sc = S2sc/100;

dataFolder = ['P' partNum '_filtered/'];
targetDist = {'55','70','85'};
i_targ = '85';
i_set = '1';
i_jump = '1';

movieName = ['P' partNum '_' i_targ '_' i_set '_' i_jump '_EKFVis_slow.avi'];

data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));


% Large Nested FOR Loops for each jump recording

%             dataName = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
%                         '_' num2str(i_jump) '_clean-P' partNum '_template'];
dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P']; %because some files names "P03-template"
file_ind = strfind(data_files,dataName);
file_ind = find(not(cellfun('isempty', file_ind)));
dataName = data_files{file_ind};

dataNameSL = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-start_land'];

%             if(exist(['../data/' dataFolder dataName '.trc'],'file')~=2)
%                 disp([dataName '.TRC does not exist'])
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
%                 trc = fillInTRCFrames(trc);


    % Read in start line and landing platform marker data
    trcSL = readTrc(['../data/' dataFolder dataNameSL '.trc']);
    markerNamesSL = fieldnames(trcSL.data);
    markerNamesSL = markerNamesSL(3:end); % first 2 names are Frame # and Time
    % Rotate markers so subject jumps in positive X direction
    for m = 1:length(markerNamesSL)
        trcSL.data.(markerNamesSL{m}) = (rotz(-pi/2)*trcSL.data.(markerNamesSL{m})')';
    end
    trcSL = fillInTRCFrames(trcSL);

    locationStart = double([mean(trcSL.data.START_L); mean(trcSL.data.START_R)])/1000; % convert to meters
    locationLand = double([mean(trcSL.data.LAND_BL); mean(trcSL.data.LAND_BR);...
                    mean(trcSL.data.LAND_TL); mean(trcSL.data.LAND_TR)])/1000;

    %Determine TO and landing frames in recording
    dt = 1/(trc.DataRate);
%                 footSpeedX = diff(trc.data.FRMed(:,1));
    footSpeedX = (diff(trc.data.FLMed(:,1)) + diff(trc.data.FRMed(:,1)))/2;
    [~,midJumpFrame] = max(footSpeedX);
    footTOFrame = find(footSpeedX(1:midJumpFrame)<=1,1,'last'); % equivalent to (1/1000)/dt = 0.2 [m/s]
    footLandFrame = midJumpFrame + find(footSpeedX(midJumpFrame:end)<=1,1,'first');
    startFrame = footTOFrame - trc.DataRate; %add 1 second to start
    endFrame = footLandFrame + trc.DataRate; %add 1 second to end
    if(endFrame > size(trc.data.FRMed,1))
        endFrame = size(trc.data.FRMed,1);
    end

%                 figure(4); clf; hold on; grid on;
%                 plot(footSpeedX,'r--');
%                 plot((footPosZ./10),'b');
%                 plot([footTOFrame,footTOFrame],[0,15],'k--');
%                 plot([footLandFrame,footLandFrame],[0,15],'k--');


    %Create Euler and Lie Group Models
    initPosWindow = [251,450];
%     [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_setLinkLengths(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath,partNum);
    [mdl_eul, trcMod_eul, initPos, modelLinks] = createJumpModel_setLinkLengths_shldrPrism(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath,partNum);

    mdl_eul.forwardPosition;
    

    if(visualize)
        vis = rlVisualizer('vis', 800, 700);
        vis.addModel(mdl_eul);
        vis.update;
    end
    
    if(makeMovie)
        vidObj = VideoWriter(movieName);
        vidObj.Quality = 100;
        vidObj.FrameRate = 20; % 200Hz data, plots every 10 frames
        open(vidObj);
    end



    disp('Position Vis camera view for movie, then continue');
    pause();
    

    %% Run EKF
    %Initialize model, set to initPos (approximate start position)
    mdl_eul.position = initPos;
    mdl_eul.velocity(:) = 0;
    mdl_eul.acceleration(:) = 0;
    mdl_eul.forwardPosition;
    mdl_eul.forwardVelocity;


    if(visualize)
        vis.update;
    end

    %Create the Measurement vector by concatenating marker data
    mes_eul = zeros(size(trc.data.Frame,1),numel(mdl_eul.sensors)*3);

    for i=1:numel(mdl_eul.sensors)
        indxs = i*3-2:i*3;
        mes_eul(:,indxs) = trcMod_eul.data.(mdl_eul.sensors(i).name);
    end


    mes_eul_predict = zeros(size(mes_eul));


    ekf_eul = EKF_Q_DQ_DDQ(mdl_eul);
    ekf_eul.observation_noise = eye(ekf_eul.sizeZ) * 0.01;
    ekf_eul.covariance = eye(ekf_eul.sizeX) * 1;
    ekf_eul.process_noise = eye(ekf_eul.sizeX) * 1.05;


    % For new marker swapping MatlabWrapper
    mes_obj = SensorMeasurement(1,30);
    [mes_obj.size] = deal(3);
    [mes_obj.type] = deal(mdl_eul.sensors(1).binary_type);
    matches = 1:numel(mdl_eul.sensors);
    matches = [matches' matches'];


    tic;
    for i=1:size(mes_eul,1)

        z_eul = mes_eul(i,:);
        mes_obj.setMesArray(z_eul); % For marker swapping MatlabWrapper
        u = dt;
        ekf_eul.run_iteration(u,mes_obj,matches); % For marker swapping MatlabWrapper


        if(visualize && (mod(i,5)==0) )
%             for m=1:numel(z_eul)/3
%                 m_pos = z_eul(m*3-2:m*3);
%                 vis.addMarker([num2str(m)],m_pos);
%             end

            vis.update();
            
            if(makeMovie)
                [img,imgbool] = vis.getScreenshot();
%                 img = flipdim(img,1);
                imshow(img);
                writeVideo(vidObj, getframe(gca)); 
            end
        end


    end
    EKF_runtime = toc



    if(visualize)
        clear vis;
    end
    
    if (makeMovie) 
        close(vidObj); 
    end

    
%     if(makeMovie)
%         mydir = pwd;
%         idcs = strfind(mydir,'\');
%         newdir = mydir(1:idcs(end)); 
%         saveFilePath = [newdir 'results\RESULTS_JA_P' partNum '_' targetDist{i_targ} '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
% 
%         world2base = mdl_eul.transforms(1).t;
% 
%         save(saveFilePath,'jointAngles_names','jointAngles_eul',...
%                 'footTOFrame','footLandFrame','startFrame','endFrame',...
%                 'modelLinks','world2base','locationStart','locationLand');
% %                     save(saveFilePath,'jointAngles_names','jointAngles_eul','jointAngles_lie',...
% %                             'footTOFrame','footLandFrame','startFrame','endFrame');
% 
%         disp([dataName ' JA data saved'])
%     end


    disp([dataName ' filtering complete'])
end


