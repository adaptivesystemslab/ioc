% Reads in previously-saved trc and ekf marker position data, replays
% visualization using ekf_state and lg_ekf_hatInvState (joint angles)

EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
addpath('..\..\');
addpath(EKFCodePath);

partNum = '04';
i_targ = '70';
i_set = '1';
i_jump = '1';

plotMarkers = 0;


data_files = cellstr(ls(fullfile(['C:\Users\kgwester\Documents\ResearchWork\Jumping Data Analysis\data\' dataFolder], '*.trc')));



mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 
loadFilePath = [newdir '\results\RESULTS_EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '.mat'];
if(exist(loadFilePath,'file')~=2)
    disp(['EKF_P' partNum '_' i_targ '_' i_set '_' i_jump '.mat NOT FOUND']);
else
    load(loadFilePath);
    % this loads in: 'mes_eul','mes_eul_predict',...
    %    'mes_lie','mes_lie_predict','ekf_state','lg_ekf_hatInvState'
    
    [gender,height,weight,S2sc,calibPose] = getPartData(partNum);
    S2sc = S2sc/100;

    % Read in .trc files and make models
    dataFolder = ['P' partNum '_filtered/'];
%     dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P' partNum '_template'];
    dataName = ['P' partNum '_target_' i_targ '_' i_set '_' i_jump '_clean-P']; %because some files names "P03-template"
    file_ind = strfind(data_files,dataName);
    file_ind = find(not(cellfun('isempty', file_ind)));
    dataName = data_files{file_ind};
    
%     dataNameSL = ['P' partNum '_target_' targetDist{i_targ} '_' num2str(i_set) ...
%                 '_' num2str(i_jump) '_clean-start_land'];

    %% Read in marker data, make models, visualize
    trc = readTrc(['../data/' dataFolder dataName]);
    markerNames = fieldnames(trc.data);
    markerNames = markerNames(3:end); % first 2 names are Frame # and Time
    % Rotate markers so subject jumps in positive X direction
    for m = 1:length(markerNames)
        trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
    end
    trc = fillInTRCFrames(trc);
    
    
    %Create Euler and Lie Group Models
    initPosWindow = [101, 200];
%             initPosWindow = [(startFrame-100), startFrame];
    [mdl_eul, trcMod_eul, initPos] = createJumpModel(trc, 1, initPosWindow, S2sc, calibPose, EKFCodePath);
    [mdl_lie, trcMod_lie, initPos_lie] = createJumpModel_Lie(trc, 1, initPosWindow, S2sc, EKFCodePath);

    mdl_eul.forwardPosition;
    %Shift the Lie group model by 1 meter in x and y to visualize better
    mdl_lie.transforms(1).t(1:2,4) = mdl_eul.transforms(1).t(1:2,4) + [1 1]';
    mdl_lie.forwardPosition();
    


    % Make Visualizer
    vis = rlVisualizer('vis',640,960);
    vis.addModel(mdl_eul);
    vis.addModel(mdl_lie);
    vis.update();
    
    %Initialize model, set to initPos (approximate start position)
    mdl_eul.position = initPos;
    mdl_eul.velocity(:) = 0;
    mdl_eul.acceleration(:) = 0;
    mdl_eul.forwardPosition;
    mdl_eul.forwardVelocity;

    mdl_lie.position = initPos_lie;
    mdl_lie.velocity(:) = 0;
    mdl_lie.acceleration(:) = 0;
    mdl_lie.forwardPosition;
    mdl_lie.forwardVelocity;

    vis.update;

    ekf = EKF_Q_DQ_DDQ(mdl_eul);
    lg_ekf = LG_EKF_Q_DQ_DDQ(mdl_lie);

    for i= 1:size(mes_eul,1)
        

        % Set EKF model state
        ekf.makeMeasure(ekf_state(i,:));
        %NOTE: "makeMeasure" function call includes model position update

        % Set LGEKF model state (getting state mapped back onto LG space)
        %FROM run_ekf_FB.m": lg_ekf.hatinvState(lg_ekf.logState(lg_ekf.state));
        lg_ekf_state = lg_ekf.expState(lg_ekf.hatState(lg_ekf_hatInvState(i,:)'));
        lg_ekf.makeMeasure(lg_ekf_state);

        if(mod(i,5)==0)
            disp(['Frame: ' num2str(i) ' out of ' num2str(size(mes_eul,1))]);
            
            if(plotMarkers)
                z_eul = mes_eul(i,:)';
                z_lie = mes_lie(i,:)';
                
                % Display Markers
                for m=1:numel(z_eul)/3
                    m_pos = z_eul(m*3-2:m*3);
                    vis.addMarker(['Meul' num2str(m)],m_pos,[1 1 0 1]);
                end
                for m=1:numel(z_lie)/3
                    m_pos = z_lie(m*3-2:m*3);
                    vis.addMarker(['Mlie' num2str(m)],m_pos,[1 1 0 1]);
                end

                
            end
            
            vis.update();
        end
        
    end
    
    clear vis;
end