classdef rlModelInstance_iit < rlModelInstance % < rlModel
    % Citations used:
    % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.872.5578&rep=rep1&type=pdf
    % https://www.tandfonline.com/doi/pdf/10.3109/17453678208992202
    
    properties
        jointCentreArray = { ...
            {'HIP_BASE_JC',         'ASIS_R', 'PSIS_R', 'ASIS_L', 'PSIS_L'}; ...
            {'SHOULDER_BASE_JC',	'SHOULDER_L', 'SHOULDER_R'}; ...
            {'HIP_R_JC',            'HARRINGTON_R'}; ...
            {'KNEE_R_JC',           'KNEE_R_LAT', 'KNEE_R_MED'}; ...
            {'ANKLE_R_JC',          'ANKLE_R_LAT', 'ANKLE_R_MED'}; ...
            {'FOOT_R_JC',           'FOOT_R_LAT', 'FOOT_R_MED'}; ...
            {'HIP_L_JC',			'HARRINGTON_L'}; ...
            {'KNEE_L_JC',			'KNEE_L_LAT', 'KNEE_L_MED'}; ...
            {'ANKLE_L_JC',			'ANKLE_L_LAT', 'ANKLE_L_MED'}; ...
            {'FOOT_L_JC',			'FOOT_L_LAT', 'FOOT_L_MED'}; ...
            {'L5S1_JC',				'ASIS_R', 'PSIS_R', 'ASIS_L', 'PSIS_L'}; ...
            {'T1C7_JC',				'SHOULDER_L', 'SHOULDER_R'}; ...0.11
            {'C1HEAD_JC',			'HEAD_T', 'HEAD_L', 'HEAD_R'}; ...
            {'SHOULDER_R_JC',		'SHOULDER_R',}; ...
            {'ELBOW_R_JC',			'ELBOW_R_LAT', 'ELBOW_R_MED',}; ...
            {'WRIST_R_JC',			'WRIST_R_LAT', 'WRIST_R_MED'}; ...
            {'HAND_R_JC',			'HANDR_T', 'HANDR_L', 'HANDR_R'}; ...
            {'SHOULDER_L_JC',		'SHOULDER_L'}; ...
            {'ELBOW_L_JC',			'ELBOW_L_LAT', 'ELBOW_L_MED'}; ...
            {'WRIST_L_JC',			'WRIST_L_LAT', 'WRIST_L_MED'}; ...
            {'HAND_L_JC',			'HANDL_T', 'HANDL_L', 'HANDL_R'}; ...
            };
        
        % link length array (linkFrame, previousFrame, previousMarker, linkMarker)
        %   linkFrame: the jointId ('fixed id' or 'revolute id') that this length will modify
        %   previousFrame: the first frame ('a idref') in the joint,
        %     coinciding with the location of previousMarker
        %   previousMarker: the marker position of previousFrame
        %   linkMarker: the second marker that describes the length
        linkLengthArray = ...
            {...
            {'length_base_rhip',                'body_base',                    'HIP_BASE_JC',     	'HIP_R_JC'}; ...	% link offsets for right lower body
            {'length_rhip_rknee',               'body_rhip_rknee',              'HIP_R_JC',         'KNEE_R_JC'}; ...
            {'length_rknee_rankle',             'body_rknee_rankle',            'KNEE_R_JC',        'ANKLE_R_JC'}; ...
            {'length_rankle_rballfoot',         'body_rankle_rballfoot',        'ANKLE_R_JC',       'FOOT_R_JC'}; ...
            {'length_base_lhip',                'body_base',                    'HIP_BASE_JC',     	'HIP_L_JC'}; ...    % link offsets for left lower body
            {'length_lhip_lknee',               'body_lhip_lknee',              'HIP_L_JC',         'KNEE_L_JC'}; ...
            {'length_lknee_lankle',             'body_lknee_lankle',            'KNEE_L_JC',        'ANKLE_L_JC'}; ...
            {'length_lankle_lballfoot',         'body_lankle_lballfoot',        'ANKLE_L_JC',       'FOOT_L_JC'}; ...
            {'length_base_l5s1',                'body_base',                    'HIP_BASE_JC',     	'L5S1_JC'}; ...     % link offsets for torso
            {'length_l5s1_t1c7',                'body_l5s1_t1c7',               'L5S1_JC',          'T1C7_JC'}; ...
            {'length_t1c7_c1head',              'body_t1c7_c1head',             'T1C7_JC',          'C1HEAD_JC'}; ...
            {'length_c7rshoulder_rshoulder',	'body_c7rshoulder_rshoulder',   'T1C7_JC',          'SHOULDER_R_JC'}; ...
            {'length_rshoulder_relbow',         'body_rshoulder_relbow',        'SHOULDER_R_JC',    'ELBOW_R_JC'}; ...
            {'length_relbow_rwrist',            'body_relbow_rwrist',           'ELBOW_R_JC',       'WRIST_R_JC'}; ...
            {'length_rwrist_rhand',            'body_rwrist_rhand',            'WRIST_R_JC',        'HAND_R_JC'}; ...
            {'length_c7lshoulder_lshoulder',    'body_c7lshoulder_lshoulder',   'T1C7_JC',          'SHOULDER_L_JC'}; ...
            {'length_lshoulder_lelbow',         'body_lshoulder_lelbow',        'SHOULDER_L_JC',    'ELBOW_L_JC'}; ...
            {'length_lelbow_lwrist',            'body_lelbow_lwrist',           'ELBOW_L_JC',       'WRIST_L_JC'}; ...
            {'length_lwrist_lhand',            'body_lwrist_lhand',            'WRIST_L_JC',       'HAND_L_JC'}; ...
            };
        
        % marker attachments (attachmentFrame, attachmentFrameMarker, [all sensorMarkers])        
        sensorAttachmentArray = ...
            {...
%             {'frame_l5s1_0',        'HIP_BASE_JC',      'ASIS_R', 'PSIS_R'}; ...      % sensor placements for right lower body
            {'frame_l5s1_0',        'HIP_BASE_JC',      'ASIS_R', 'PSIS_R', 'ASIS_L', 'PSIS_L'}; ...      % sensor placements for right lower body
            {'frame_rknee_0',       'KNEE_R_JC',        'KNEE_R_LAT', 'KNEE_R_MED'}; ...               
            {'frame_rankle_0',      'ANKLE_R_JC',       'ANKLE_R_LAT', 'ANKLE_R_MED'}; ...
            {'frame_rballfoot_end', 'FOOT_R_JC',        'FOOT_R_LAT', 'FOOT_R_MED', 'HEEL_R'}; ... % 
            {'frame_lknee_0',       'KNEE_L_JC',        'KNEE_L_LAT', 'KNEE_L_MED'}; ...                 % sensor placements for left lower body
            {'frame_lankle_0',      'ANKLE_L_JC',       'ANKLE_L_LAT', 'ANKLE_L_MED'}; ...
            {'frame_lballfoot_end', 'FOOT_L_JC',        'FOOT_L_LAT', 'FOOT_L_MED', 'HEEL_L'}; % 
%             ...
%             {'frame_t1c7_0',        'T1C7_JC',          'C7'};
            {'frame_rshoulder_0',   'SHOULDER_R_JC',    'SHOULDER_R'}; ...                              % sensor placements for right upper body
            {'frame_relbow_0',      'ELBOW_R_JC',       'ELBOW_R_LAT', 'ELBOW_R_MED'}; ...
            {'frame_rwrist_0',      'WRIST_R_JC',       'WRIST_R_LAT', 'WRIST_R_MED'}; ...
            {'frame_lshoulder_0',   'SHOULDER_L_JC',    'SHOULDER_L'}; ...                              % sensor placements for left upper body
            {'frame_lelbow_0',      'ELBOW_L_JC',       'ELBOW_L_LAT', 'ELBOW_L_MED'}; ...
            {'frame_lwrist_0',      'WRIST_L_JC',       'WRIST_L_LAT','WRIST_L_MED'}; ...
            };
        
        % trc/imu attachments (attachmentFrame, [name for the cluster assocation], [all sensorMarkers])
        sensorSecondaryAttachmentArray = ...
            {...
            {'frame_l5s1_0',        'HIP_BASE_JC',      'PELV_T', 'PELV_L', 'PELV_R'}; ...
            {'frame_rknee_0',       'KNEE_R_JC',       	'uLEGR_T', 'uLEGR_L', 'uLEGR_R'}; ...
            {'frame_rankle_0',      'ANKLE_R_JC',       'lLEGR_T', 'lLEGR_L', 'lLEGR_R'}; ...
            {'frame_rballfoot_end', 'FOOT_R_JC',        'FOOTR_T',  'FOOTR_L', 'FOOTR_R'}; ...
            {'frame_lknee_0',       'KNEE_L_JC',        'uLEGL_T', 'uLEGL_L', 'uLEGL_R'}; ...
            {'frame_lankle_0',      'ANKLE_L_JC',       'lLEGL_T', 'lLEGL_L', 'lLEGL_R'}; ...
            {'frame_lballfoot_end', 'FOOT_L_JC',        'FOOTL_T', 'FOOTL_L', 'FOOTL_R'}; ...
% %             ...
            {'frame_t1c7_0',        'T1C7_JC',          'STERN_T', 'STERN_L', 'STERN_R'}; ...           % sensor placement for all the XSEN IMUs
            {'frame_t1c7_0',        'T1C7_JC',          'HEAD_T', 'HEAD_L', 'HEAD_R'}; ...
            {'frame_t1c7_0',        'T1C7_JC',          'SHOUR_T', 'SHOUR_L', 'SHOUR_R'}; ...
            {'frame_relbow_0',      'ELBOW_R_JC',       'uARMR_T', 'uARMR_L', 'uARMR_R'}; ...
            {'frame_rwrist_0',      'WRIST_R_JC',       'fARMR_T', 'fARMR_L', 'fARMR_R'}; ...
            {'body_rwrist_rhand',   'WRIST_R_JC',       'HANDR_T', 'HANDR_L', 'HANDR_R'}; ...
            {'frame_t1c7_0',        'T1C7_JC',          'SHOUL_T', 'SHOUL_L', 'SHOUL_R'}; ...
            {'frame_lelbow_0',      'ELBOW_L_JC',       'uARML_T', 'uARML_L', 'uARML_R'}; ...
            {'frame_lwrist_0',      'WRIST_L_JC',       'fARML_T', 'fARML_L', 'fARML_R'}; ...
            {'body_lwrist_lhand',   'WRIST_L_JC',       'HANDL_T', 'HANDL_L', 'HANDL_R'}; ...
            };       
        
		% joint limits [deg]
        jointRangeOfMotion = {};
%         jointRangeOfMotion = {...
% %             {'joint_rshoulder_0',	[-62,	166]}; ...
% %             {'joint_rshoulder_1',	[-184,	0]}; ...
% %             {'joint_rshoulder_2', 	[-68,	103]}; ...
%             {'joint_relbow_0',		[0,		142]   }; ...
% %             {'joint_relbow_1',		[-82+90,75+90]}; ...
%             {'joint_rwrist_0',		[-75,  	76]   }; ...
% %             {'joint_rwrist_1',		[-36,  	21]   }; ...
%             ...
% %             {'joint_rhip_0',		[-9, 	122]}; ...
% %             {'joint_rhip_1',		[-45,  	26]}; ...
% %             {'joint_rhip_2',		[-47,  	47]}; ...
%             {'joint_rknee_0',		[-142,  1]}; ...
%             {'joint_rankle_0',		[-56, 	12]}; ...
% %             {'joint_rankle_2',      [-18, 	33]}; ...
%             ...
% %             {'joint_lshoulder_0',  	[-62, 	166]}; ...
% %             {'joint_lshoulder_1',  	[0, 	184]}; ...
% %             {'joint_lshoulder_2', 	[103+90,68+270]}; ...
%             {'joint_lelbow_0',      [0, 	142]}; ...
% %             {'joint_lelbow_1',     	[75+90, 82+270]}; ...
%             {'joint_lwrist_0',    	[-75,  	76]   }; ...
% %             {'joint_lwrist_1',    	[-21,  	36]   }; ...
%             ...
% %             {'joint_lhip_0',    	[-9, 	122]}; ...
% %             {'joint_lhip_1',    	[-26,  	45]}; ...
% %             {'joint_lhip_2',     	[-47,  	47]}; ...
%             {'joint_lknee_0',    	[-142,  1]}; ...
%             {'joint_lankle_0',      [-56, 	12]}; ...
% %             {'joint_lankle_2',		[-33, 	18]}; ...
%             };
        
        jointRangeOfMotionTolerance = 30;
        
        % subject demographics (id, height [m], weight [kg])
        subjectDemographicArray = {...
            {1, 1.70, 62, 'f'} ...
            {2, 1.70, 70, 'm'} ...
			{3, 1.69, 62, 'm'} ...
			{4, 1.71, 57, 'm'} ...
			{5, 1.62, 60, 'f'} ... % uncertain weight
			{6, 1.84, 75, 'm'} ... % uncertain weight
			{7, 1.76, 76, 'm'} ...
			{8, 1.80, 77, 'm'} ...
			{9, 1.78, 75, 'm'} ...
			{10, 1.82, 78, 'm'} ... % uncertain weight
			{11, 1.65, 64, 'f'} ...
			{12, 1.79, 83, 'm'} ...
			{13, 1.74, 80, 'm'} ... % uncertain weight
			{14, 1.86, 94, 'm'} ...
			{15, 1.87, 99, 'm'} ...
            };

        linkLengthSymmetry = {...
            {'length_base_rhip', 'length_base_lhip'}, ...
            {'length_rhip_rknee', 'length_lhip_lknee'}, ...
            {'length_rknee_rankle', 'length_lknee_lankle'}, ...
            {'length_rankle_rballfoot', 'length_lankle_lballfoot'}, ...
            {'length_c7rshoulder_rshoulder', 'length_c7lshoulder_lshoulder'}, ...
            {'length_rshoulder_relbow', 'length_lshoulder_lelbow'}, ...
            {'length_relbow_rwrist', 'length_lelbow_lwrist'}, ...
            {'length_rwrist_rhand', 'length_lwrist_lhand'}, ...
            };
        
        sensorSymmetry = {...
            {'KNEE_R_LAT', 'KNEE_R_MED', 'KNEE_L_LAT', 'KNEE_L_MED'}, ...
            {'ANKLE_R_LAT', 'ANKLE_R_MED', 'ANKLE_L_LAT', 'ANKLE_L_MED'}, ...
            {'FOOT_R_LAT', 'FOOT_R_MED', 'FOOT_L_LAT', 'FOOT_L_MED'}, ...
            {'HEEL_R', 'HEEL_L'}, ...
            {'SHOULDER_R', 'SHOULDER_L'}, ...
            {'ELBOW_R_LAT', 'ELBOW_R_MED', 'ELBOW_L_LAT', 'ELBOW_L_MED'}, ...
            {'WRIST_R_LAT', 'WRIST_R_MED', 'WRIST_L_LAT','WRIST_L_MED'}, ...
            };

        lengthUseInd = 1:10;

        harrington_crop1 = 1;
        harrington_crop2 = 10;
        
        trcSensorDecorator = {'position'};
        imuSensorDecorator = {'gyroscope', 'accelerometer'};
        imuSensorYaw = 1;

        mocapMarkers_rightArm = {'SHOULDER_R', 'ELBOW_R_LAT', 'ELBOW_R_MED', 'WRIST_R_LAT', 'WRIST_R_MED'};
        mocapMarkers_leftArm = {'SHOULDER_L', 'ELBOW_L_LAT', 'ELBOW_L_MED', 'WRIST_L_LAT', 'WRIST_L_MED'};
        mocapMarkers_rightLeg = {'KNEE_R_LAT', 'KNEE_R_MED', 'ANKLE_R_LAT', 'ANKLE_R_MED', 'FOOT_R_LAT', 'FOOT_R_MED', 'HEEL_R'};
        mocapMarkers_leftLeg = {'KNEE_L_LAT', 'KNEE_L_MED', 'ANKLE_L_LAT', 'ANKLE_L_MED', 'FOOT_L_LAT', 'FOOT_L_MED', 'HEEL_L'};
        mocapMarkers_torso = {'ASIS_R', 'PSIS_R', 'ASIS_L', 'PSIS_L', 'C7'};
        
        joints_rightArm = {'joint_c7rshoulder_0', 'joint_rshoulder_0', 'joint_rshoulder_1', 'joint_rshoulder_2', 'joint_relbow_0', 'joint_relbow_1'};
        joints_leftArm = {'joint_c7lshoulder_0', 'joint_lshoulder_0', 'joint_lshoulder_1', 'joint_lshoulder_2', 'joint_lelbow_0', 'joint_lelbow_1'};
        joints_rightLeg = {'joint_rhip_0', 'joint_rhip_1', 'joint_rhip_2', 'joint_rknee_0', 'joint_rankle_0', 'joint_rankle_1', 'joint_rankle_2'};
        joints_leftLeg = {'joint_lhip_0', 'joint_lhip_1', 'joint_lhip_2', 'joint_lknee_0', 'joint_lankle_0', 'joint_lankle_1', 'joint_lankle_2'};
        joints_torso = {'frame_6dof_revz0_to_frame_6dof_revz1', 'frame_6dof_revy0_to_frame_6dof_revy1', 'frame_6dof_revx0_to_frame_6dof_revx1'};
        
        jointsPlot = {'joint_rhip_0', 'joint_rhip_1', 'joint_rhip_2', 'joint_rknee_0', ...
            'joint_lhip_0', 'joint_lhip_1', 'joint_lhip_2', 'joint_lknee_0', ...
            'frame_6dof_revz0_to_frame_6dof_revz1', 'frame_6dof_revy0_to_frame_6dof_revy1', 'frame_6dof_revx0_to_frame_6dof_revx1'};
    end
    
    methods
        function obj = rlModelInstance_iit(subjectId)
            obj.datasetName = 'IIT';
            obj.subjectId = subjectId;
            obj.loadDemographics();
        end
        
        function [dt, time, data] = loadData_trc(obj, pathToTrc, algorithmParam)
            % check to see if there is any 'trim' version of the data. if there
            % is, use that one instead, and iterate though all the trims
            [pathstr,name,ext] = fileparts(pathToTrc); 
            trcPrefix = fullfile(pathstr, [name '_trim*.trc']);
            dirTrc = dir(trcPrefix);
            
            pathsToLoad = {};
            if isempty(dirTrc)
                pathsToLoad{1} = pathToTrc;
            else
                for i = 1:length(dirTrc)
                    currFile = dirTrc(i).name;
                    if length(strsplit(currFile, 'Unnamed')) > 1
                        continue
                    end
                    
                    pathsToLoad{end+1} = fullfile(pathstr, currFile);
                end
            end
        
            numUnknownMarkersLoaded = 0;
            for i = 1:length(pathsToLoad)
%                 if i > 1
%                     for j = 1:10
%                         fprintf('rlModelInstance_iit: WARNING current configuration cuts off at 1 trims\n');
%                     end
%                     continue;
%                 end
                
                % load the base file
                trcPath = pathsToLoad{i};
                trc_named = loadDataFromTrc_iit(trcPath, algorithmParam);
                trc_named = fillInTRCFrames(trc_named, obj.lengthUseInd(end));

                % also load the unlabelled markers into the file
                [pathstr,name,ext] = fileparts(trcPath);
                trcPath = fullfile(pathstr, [name '-Unnamed' ext]);
                trc_unnamed = loadDataFromTrc_iit(trcPath, algorithmParam);
                
                % merge the unnamed markers into the named ones
                field_names = fieldnames(trc_unnamed.data);
                for j = 1:numel(field_names)
                    % concat current entries into the previous one
                    numUnknownMarkersLoaded = numUnknownMarkersLoaded + 1;
                    sensorName = field_names{j};
                    
                    % note that we are not assuming that unknown markers over
                    % different mocap trims are continous. to assume they are
                    % continous, just set newSensorName to sensorName
%                     newSensorName = ['U' pad(num2str(numUnknownMarkersLoaded), 8, 'left', '0')];
                    newSensorName = sensorName;
                    trc_named.data.(newSensorName) = trc_unnamed.data.(sensorName);
                end
                
                if i == 1
                    trc = trc_named;
                else % concat new trcdata entries into the previous one
                    indStart = length(trc.time) + 1;
                    indEnd = indStart + length(trc_named.time) - 1;
                    inds = indStart:indEnd;
                    
                    % start with the time array
                    newTime = trc.time(end) + trc.dt + trc_named.time;
                    trc.time(inds, :) = newTime;

                    % prepend each trcdata with zeros first. this way, if the
                    % new trcnamed does not have (sensorname) that trcdata has,
                    % the array lengths will not become misaligned
                    field_names = fieldnames(trc.data); 
                    newArray = zeros(size(inds, 2), 3);
                    for j = 1:numel(field_names)
                        sensorName = field_names{j};
                        trc.data.(sensorName)(inds, :) = newArray;
                    end
                    
                    % then add the data
                    field_names = fieldnames(trc_named.data); 
                    for j = 1:numel(field_names)
                        sensorName = field_names{j};
      
                        % if the original trcdata does not have (sensorName),
                        % prepend one for it and fill it with zeros
                        if ~isfield(trc.data, sensorName)
                            trc.data.(sensorName) = zeros(size(trc.time, 1), 3);
                        end
                        
                        trc.data.(sensorName)(inds, :) = trc_named.data.(sensorName);
                    end
                end
            end
            
            dt = mean(diff(trc.time));
            time = trc.time;
            data = trc.data;
        end
        
        function [kinematicTransform, dynamicTransform] = applyLinkTransform(obj, targetFrameStr, sourceFrameStr, linkVectorMean, lengthUseInd, dataInstance)
            % variable setup
            rotMirror = [-1 0 0; 0 1 0; 0 0 1];

            linkLength = norm(linkVectorMean);
            
            % set any link-specific information
            switch targetFrameStr
                case {'length_base_l5s1', 'length_base_lhip', ...
                        'length_c7rshoulder_rshoulder', 'length_c7lshoulder_lshoulder'} % either doesn't exist or already accounted for in other limbs
                    dumasFrameStr = '';
                    
                case {'length_base_rhip'}
                    dumasFrameStr = 'pelvis';
                    
                case {'length_rhip_rknee', 'length_lhip_lknee'}
                    dumasFrameStr = 'thigh';
                    
                case {'length_rknee_rankle', 'length_lknee_lankle'}
                    dumasFrameStr = 'leg';
                    
                case {'length_rankle_rballfoot', 'length_lankle_lballfoot'}
                    dumasFrameStr = 'foot';
                    
                case {'length_l5s1_t1c7'}
                    dumasFrameStr = 'torso';
                    
                case {'length_t1c7_c1head'}
                    dumasFrameStr = 'head&neck';
                    
                case {'length_rshoulder_relbow', 'length_lshoulder_lelbow'}
                    dumasFrameStr = 'arm';
                    
                case {'length_relbow_rwrist', 'length_lelbow_lwrist'}
                    dumasFrameStr = 'forearm';
                    
                case {'length_rwrist_rhand', 'length_lwrist_lhand'}
                    dumasFrameStr = 'hand';
            end
            
            if ~isempty(dumasFrameStr)
                anthLinkDistance =            lookupTableDumas('length',     dumasFrameStr, obj.gender, [], [])*obj.body_height;
            else
                anthLinkDistance = 0;
            end

            switch targetFrameStr 
                case {'length_base_rhip'}
                    X00Length = linkVectorMean;
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2);  % apply the full pelvis to the right side instead of splitting it up
                    anthLength = [anthLinkDistance 0 0]';

                case {'length_base_lhip'}
                    X00Length = linkVectorMean;
                    rotMtxDumas = []; 
                    anthLength = [-anthLinkDistance 0 0]';

                case {'length_c7rshoulder_rshoulder'}
                    X00Length = [linkLength 0 0]';
                    rotMtxDumas = [];
%                     [X00Length2, ~] = amabile2006centre(dataInstance.data, obj.harrington_crop1, obj.harrington_crop2, obj.body_height);
           
                case {'length_c7lshoulder_lshoulder'}
                    X00Length = [-linkLength 0 0]';
                    rotMtxDumas = [];
                    
                case {'length_rankle_rballfoot'}
                    X00Length = [0 linkLength 0]';
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2);
                    anthLength = [0 anthLinkDistance 0]';

                 case {'length_lankle_lballfoot'}
                    X00Length = [0 linkLength 0]';
                    rotMtxDumas = rotMirror*rotz(pi/2)*rotx(pi/2); 
                    anthLength = [0 anthLinkDistance 0]';

                case {'length_base_l5s1', 'length_l5s1_t1c7'} % first joint is sagittal
                    X00Length = [0 0 linkLength]';
                    rotMtxDumas = roty(pi)*rotz(pi/2)*rotx(pi/2);  % rotating from Dumas, but pointing up
                    anthLength = [0 0 anthLinkDistance]';

                case {'length_t1c7_c1head'}
%                     X00Length = [0 0 linkLength]';
                    X00Length = [0 0 0]';
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2);
                    anthLength = [0 0 0]';

                case {'length_rhip_rknee', 'length_rknee_rankle'} % first joint is sagittal
                    X00Length = [0 0 -linkLength]';
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2); % dumas down
                    anthLength = [0 0 -anthLinkDistance]';

                case {'length_lhip_lknee', 'length_lknee_lankle'} % first joint is sagittal
                    X00Length = [0 0 -linkLength]';
                    rotMtxDumas = rotMirror*rotz(pi/2)*rotx(pi/2); % dumas down
                    anthLength = [0 0 -anthLinkDistance]';
                
                case {'length_rshoulder_relbow', 'length_relbow_rwrist', 'length_rwrist_rhand'}
                    X00Length = [0 0 -linkLength]';
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2); % dumas down
                    anthLength = [0 0 -anthLinkDistance]';

                case {'length_lshoulder_lelbow', 'length_lelbow_lwrist', 'length_lwrist_lhand'}
                    X00Length = [0 0 -linkLength]';
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2); % dumas down
                    anthLength = [0 0 -anthLinkDistance]';
            end
                  
            switch obj.linkDefinition
                case 'initialpose'
                    correctedLinkLength = linkLength;
                 
                case 'X00'
                    correctedLinkLength = X00Length;
                    
%                     if isnan(linkLength)
%                         correctedLinkLength = anthLength;
%                     end
            end

            T = eye(4);
            T(1:3, 4) = correctedLinkLength;
          
            if ~isempty(dumasFrameStr)
                mass =            lookupTableDumas('mass',     dumasFrameStr, obj.gender, [], [])*obj.body_weight;
                comScale =        lookupTableDumas('com',      dumasFrameStr, obj.gender, linkLength, []);
                inertialScale =   lookupTableDumas('inertial', dumasFrameStr, obj.gender, linkLength, mass);
                
                % parallel axis, since Dumas inertia assumes it's on the edge of
                % the limb
                com = rotMtxDumas*comScale;
                inertial = rotMtxDumas*inertialScale*rotMtxDumas';
%                 inertial = zeros(3, 3);
                   
%                 targetFrameStr
%                 [maxVal, maxInd] = max(abs([linkLength X00Length com]));
%                 vec = [linkLength X00Length comScale com];
            else
                mass = [];
                com = [];
                inertial = [];
            end
 
            [kinematicTransform, dynamicTransform] = obj.assembleKinDynTransform(targetFrameStr, sourceFrameStr, T, mass, com, inertial);
        end
        
        function sensorTransform = applyMarkerSensor(obj, sourceFrameStr, targetMarkerStr, linkVectorMean, lengthUseInd, dataInstance)   
            linkDistance = norm(linkVectorMean);
            
            switch targetMarkerStr
% % %                 case {'ASIS_R', 'ASIS_L', 'PSIS_L', 'PSIS_R'} 
% % %                     switch sensorMarkerStr
% % %                         case {'ASIS_R', 'ASIS_L'} 
% % %                               linkLengthArrangement_R = {sourceFrameStr, 'HIP_BASE_JC', 'ASIS_R'};
% % %                               linkLengthArrangement_L = {sourceFrameStr, 'HIP_BASE_JC', 'ASIS_L'};
% % %                               
% % %                         case {'PSIS_L', 'PSIS_R'} 
% % %                              linkLengthArrangement_R = {sourceFrameStr, 'HIP_BASE_JC', 'PSIS_R'};
% % %                              linkLengthArrangement_L = {sourceFrameStr, 'HIP_BASE_JC', 'PSIS_L'};
% % %                     end
% % % 
% % %                     [attachmentFrameStr_test, attachmentMarkerStr_test, sensorMarkerStr_test, markerOffset_R, sensorMarker_test, sensorMarkerMean_test, sensorMarkerStd_test] = ...
% % %                         obj.sensorAttachmentArrayToMarkerData(linkLengthArrangement_R, lengthUseInd, dataInstance);
% % %                     checkSensBool_R = checkSensorPlacement([], sensorMarker_test, sensorMarkerMean_test, sensorMarkerStd_test, ...
% % %                         markerOffset_R, lengthUseInd, sensorMarkerStr_test, attachmentFrameStr_test);
% % %                     
% % %                     [attachmentFrameStr_test, attachmentMarkerStr_test, sensorMarkerStr_test, markerOffset_L, sensorMarker_test, sensorMarkerMean_test, sensorMarkerStd_test] = ...
% % %                         rlModelInstance.sensorAttachmentArrayToMarkerData(linkLengthArrangement_L, lengthUseInd, dataInstance);
% % %                     checkSensBool_L = checkSensorPlacement([], sensorMarker_test, sensorMarkerMean_test, sensorMarkerStd_test, ...
% % %                         markerOffset_L, lengthUseInd, sensorMarkerStr_test, attachmentFrameStr_test);
% % %                         
% % %                     markerOffsetSign = sign(markerOffset);
% % %                     meanAbsMagitude = [markerOffset_R; markerOffset_L];
% % %                     if checkSensBool_R == 1 && checkSensBool_R == 0
% % %                         meanAbsMagitude = [markerOffset_R];
% % %                     elseif checkSensBool_R == 0 && checkSensBool_R == 1
% % %                         meanAbsMagitude = [markerOffset_L];
% % %                     end
% % % 
% % %                     X00Length = markerOffsetSign .* mean(abs(meanAbsMagitude), 1); % ensure hip markers are symmetric
  
                case 'C7'
                    X00Length = [0 -linkDistance 0]';

                case {'SHOULDER_R', 'SHOULDER_L'} % no offset
                    X00Length = [0 0 0]';
                  
                case {'ELBOW_R_MED', 'WRIST_R_MED', ...
                        'ELBOW_L_LAT', 'WRIST_L_LAT'} % left of the link
                    X00Length = [-linkDistance 0 0]';

                case {'ELBOW_R_LAT', 'WRIST_R_LAT', ...
                        'ELBOW_L_MED', 'WRIST_L_MED'} % right of the link
                    X00Length = [linkDistance 0 0]';

                case {'KNEE_R_MED', 'ANKLE_R_MED', 'FOOT_R_MED', ...
                        'KNEE_L_LAT', 'ANKLE_L_LAT', 'FOOT_L_LAT'} % left of the link
                    X00Length = [-linkDistance 0 0]';

                case {'KNEE_R_LAT', 'ANKLE_R_LAT', 'FOOT_R_LAT', ...
                        'KNEE_L_MED', 'ANKLE_L_MED', 'FOOT_L_MED'} % right of the link
                    X00Length = [linkDistance 0 0]';
                    
                case {'HEEL_R', 'HEEL_L'} % right of the link
                    X00Length = [0 -linkDistance 0]';

                otherwise % IMU markers
                    X00Length = linkVectorMean;
            end
            
            switch obj.linkDefinition
                case 'initialpose'
                    correctedLinkLength = linkVectorMean;
                    
                case 'X00'
                    correctedLinkLength = X00Length;
            end
            
            T = eye(4);
            T(1:3,4) = correctedLinkLength;
            
            sensorTransform = obj.assembleSenTransform(targetMarkerStr, sourceFrameStr, obj.trcSensorDecorator, T);
        end
        
        function obj = loadData_fp(obj, pathtoAnc, pathToAncUnloaded, algorithmParam)
            obj.fpData = loadDataFromAnc_iit(pathtoAnc, pathToAncUnloaded, algorithmParam);
        end
    end
end
