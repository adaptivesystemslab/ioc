classdef rlModelInstance_expressiveioc_rightArm < rlModelInstance % < rlModel
    properties
%         jointCentreArray = { ...
%             {'SHOULDER_R_JC',		'RCAJ',}; ...
%             {'ELBOW_R_JC',			'RHLE', 'RHME',}; ...
%             {'WRIST_R_JC',			'RUSP', 'RRSP'}; ...
%             {'HAND_R_JC',			'RHA1',  'RHA2'}; ...
%             };
				
        jointCentreArray = { ...         
            {'SHOULDER_R_JC',		'RCAJ',}; ...
            {'ELBOW_R_JC',			'RHLE', 'RHME',}; ...
            {'WRIST_R_JC',			'RUSP', 'RRSP'}; ...
            {'HAND_R_JC',			'RHA1',  'RHA2', 'RHA3', 'RHA4'}; ...
            };
        
        % link length array (linkFrame, previousFrame, previousMarker, linkMarker)
        %   linkFrame: the jointId ('fixed id' or 'revolute id') that this length will modify
        %   previousFrame: the first frame ('a idref') in the joint,
        %     coinciding with the location of previousMarker
        %   previousMarker: the marker position of previousFrame
        %   linkMarker: the second marker that describes the length
        linkLengthArray = ...
            {...
            {'length_rshoulder_relbow',         'body_rshoulder_relbow',        'SHOULDER_R_JC',    'ELBOW_R_JC'}; ...
            {'length_relbow_rwrist',            'body_relbow_rwrist',           'ELBOW_R_JC',       'WRIST_R_JC'}; ...
            {'length_rwrist_rhand',             'body_rwrist_rhand',            'WRIST_R_JC',       'HAND_R_JC'}; ...
            };

        % marker attachments (attachmentFrame, attachmentFrameMarker, [all sensorMarkers])        
        % Joint definition markers
        sensorAttachmentArray = ...
            {...
            {'frame_rshoulder_0',   'SHOULDER_R_JC',    'RCAJ'}; ...                              % sensor placements for right upper body
            {'frame_relbow_0',      'ELBOW_R_JC',       'RHLE', 'RHME'}; ...
            {'frame_rwrist_0',      'WRIST_R_JC',       'RUSP', 'RRSP'}; ...
            };
        
        % trc/imu attachments (attachmentFrame, [name for the cluster assocation], [all sensorMarkers])
        % Tracking markers
%         sensorSecondaryAttachmentArray = ...
%             {...
%             {'frame_rshoulder_6',   'SHOULDER_R_JC',    'RUPA1', 'RUPA2', 'RUPA3', 'RUPA4'};...
%             {'frame_relbow_4',      'ELBOW_R_JC',       'RFA1', 'RFA2', 'RFA3', 'RFA4'}; ...
%             {'frame_rwrist_4',      'WRIST_R_JC',       'RHA1', 'RHA2', 'RHA3', 'RHA4'}; ...
%             %{'frame_6dof_root',     '',                 'LAH', 'LPH', 'RPH', 'RAH'}; ...
%             };       
%         sensorSecondaryAttachmentArray = ...
%             {...
%             {'frame_l5_0',        'L5_JC',      'LIAS', 'RIAS', 'RIPS', 'LIPS'}; ...
%             {'frame_l5_0',        'L5_JC',      'LUS1', 'LUS2', 'LUS3', 'LUS4'}; ...
%             {'frame_t10_0',        'T10_JC',      'CV7', 'TV10'}; ...
%             {'frame_t10_0',        'T10_JC',          'THS1', 'THS2', 'THS3', 'THS4'}; ...           % sensor placement for all the XSEN IMUs
%             {'frame_c7_0',        'C7_JC',          'LAH', 'LPH', 'RPH', 'RAH'}; ...
%             {'frame_c7_0',        'C7_JC',          'SJN', 'SXS'}; ... 
%             {'frame_rshoulder_0',   'SHOULDER_R_JC',    'RUPA1', 'RUPA2', 'RUPA3', 'RUPA4'};...
%             {'frame_relbow_0',      'ELBOW_R_JC',       'RFA1', 'RFA2', 'RFA3', 'RFA4'}; ...
%             {'frame_rwrist_0',      'WRIST_R_JC',       'RHA1', 'RHA2', 'RHA3', 'RHA4'}; ...
%             {'frame_lshoulder_0',   'SHOULDER_L_JC',    'LUPA4', 'LUPA3', 'LUPA2', 'LUPA1'};...
%             {'frame_lelbow_0',      'ELBOW_L_JC',       'LFA4', 'LFA3', 'LFA2', 'LFA1'}; ...
%             {'frame_lwrist_0',      'WRIST_L_JC',       'LHA4', 'LHA3', 'LHA2', 'LHA1'}; ...
%             };       
        
        sensorSecondaryAttachmentArray = ...
            {...
            {'frame_relbow_0',      'ELBOW_R_JC',   'RUPA1', 'RUPA2', 'RUPA3', 'RUPA4'};...
            {'frame_rwrist_0',      'WRIST_R_JC',   'RFA1', 'RFA2', 'RFA3', 'RFA4'}; ...
            {'frame_rhand_end',  	'HAND_R_JC',    'RHA1', 'RHA2', 'RHA3', 'RHA4'}; ...
            };

		% joint limits [deg]
        jointRangeOfMotion = {};
        jointRangeOfMotionTolerance = 30;
        
        % subject demographics (id, height [m], weight [kg])
        subjectDemographicArray = {...
            {1, 1.58, 57, 'f'} ...
        };

        linkLengthSymmetry = {};
        sensorSymmetry = {};

        lengthUseInd = 1:10;

        harrington_crop1 = 1;
        harrington_crop2 = 10;
        
        trcSensorDecorator = {'position'};
    end
    
    methods
        function obj = rlModelInstance_expressiveioc_rightArm(subjectId)
            obj.datasetName = 'IOC';
            obj.subjectId = subjectId;
            obj.loadDemographics();
        end
        
        function [dt, time, data] = loadData_trc(obj, pathToTrc, algorithmParam)
            % check to see if there is any 'trim' version of the data. if there
            % is, use that one instead, and iterate though all the trims 
            trc = loadDataFromTrc_expressiveioc(pathToTrc, algorithmParam);
            trc = fillInTRCFrames(trc, obj.lengthUseInd(end));
            
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
                case {'length_base_l5', 'length_base_lhip', 'length_l5_t10', ...
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

                case {'length_t10_c7'}
                    dumasFrameStr = 'torso';
                    
                case {'length_c7_c1head'}
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
%                case {'length_rshoulder_relbow', 'length_relbow_rwrist', 'length_rwrist_rhand'}
%                    X00Length = [linkLength 0 0]';    
%                     X00Length = [1 0 0]';
%                    rotMtxDumas = rotz(pi/2); % dumas down

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

                case {'length_base_l5', 'length_l5_t10', 'length_t10_c7'} % first joint is sagittal
                    X00Length = [0 0 linkLength]';
                    rotMtxDumas = roty(pi)*rotz(pi/2)*rotx(pi/2);  % rotating from Dumas, but pointing up
                    anthLength = [0 0 anthLinkDistance]';

                case {'length_c7_c1head'}
%                     X00Length = [0 0 linkLength]';
                    X00Length = [0 0 linkLength]';
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2);
                    anthLength = [0 0 0]';

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
                    correctedLinkLength = linkVectorMean;
                 
                case 'X00'
                    correctedLinkLength = X00Length;
                    
                    if isnan(linkLength)
                        correctedLinkLength = anthLength;
                    end
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
                case {'RCAJ', 'LCAJ'} % no offset
                    X00Length = [0 0 0]';
                    
%                 case {'CV7', 'SJN', 'SXS', 'TV10'}
%                     X00Length = [0 linkVectorMean(2) linkVectorMean(3)];
%                   
                case {'RHLE', 'RUSP', 'LHME', 'LRSP'} 
                    X00Length = [linkDistance 0 0]';

                case {'RHME', 'RRSP', 'LHLE', 'LUSP'} 
                    X00Length = [-linkDistance 0 0]';

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