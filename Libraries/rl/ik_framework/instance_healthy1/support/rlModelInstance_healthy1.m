classdef rlModelInstance_healthy1 < rlModelInstance % < rlModel
    properties
        jointCentreArray = { ...
%             {'SHOULDER_R_JC',	'SHOULDER_R'}; ...
            {'HIP_R_JC',        'ASIS_R'}; ...
            {'KNEE_R_JC',       'KNEE_R_LAT', 'KNEE_R_MED'}; ...
            {'ANKLE_R_JC',      'ANKLE_R_LAT', 'ANKLE_R_MED'}; ...     
        };
        
        % link length array (linkFrame, previousFrame, previousMarker, linkMarker)
        %   linkFrame: the jointId ('fixed id' or 'revolute id') that this length will modify
        %   previousFrame: the first frame ('a idref') in the joint,
        %     coinciding with the location of previousMarker
        %   previousMarker: the marker position of previousFrame
        %   linkMarker: the second marker that describes the length
        linkLengthArray = ...
            {...
%             {'length_rshoulder_rhip',   'body_rshoulder_rhip',   'SHOULDER_R_JC',  'HIP_R_JC'}; ...
            {'length_rhip_rknee',   'body_rhip_rknee',   'HIP_R_JC',  'KNEE_R_JC'}; ...
            {'length_rknee_rankle', 'body_rknee_rankle', 'KNEE_R_JC', 'ANKLE_R_JC'}; ...
            };
        
        % marker arrangements (attachmentFrame, attachmentFrameMarker, [all sensorMarkers])
%         sensorAttachmentArray = ...
%             {...
%             {'frame_rshoulder_0',    'SHOULDER_R_JC', 'SHOULDER_R'}; ...
%             {'frame_rhip_0',    'HIP_R_JC',   'ASIS_R'}; ...
%             {'frame_rknee_0',   'KNEE_R_JC',  'KNEE_R_MED', 'KNEE_R_LAT'}; ...
%             {'frame_rankle_0',   'ANKLE_R_JC', 'ANKLE_R_MED', 'ANKLE_R_LAT'}; ...
%             };
        sensorAttachmentArray = ...
            {...
%             {'frame_rshoulder_0',  	'SHOULDER_R_JC',	'SHOULDER_R', 'SHOULDER_L'}; ...
            {'frame_rhip_0',        'HIP_R_JC',         'ASIS_R', 'ASIS_L'}; ...
            {'frame_rknee_0',       'KNEE_R_JC',        'KNEE_R_MED', 'KNEE_R_LAT'}; ...
            {'frame_rankle_0',      'ANKLE_R_JC',       'ANKLE_R_MED', 'ANKLE_R_LAT'}; ...
%             {'frame_rknee_0',   'KNEE_R_JC',  'KNEE_R_JC', 'KNEE_R_IMU'}; ...
%             {'frame_rankle_0',   'ANKLE_R_JC', 'ANKLE_R_JC', 'ANKLE_R_IMU'}; ...
%             {'frame_rknee_0',   'KNEE_R_JC',  'KNEE_R_MED', 'KNEE_R_LAT', 'KNEE_R_IMU'}; ...
%             {'frame_rankle_0',   'ANKLE_R_JC', 'ANKLE_R_MED', 'ANKLE_R_LAT', 'ANKLE_R_IMU'}; ...
            };
        
        sensorSecondaryAttachmentArray = ...
            {...
%             {'body_rshoulder_rhip', 'HIP_R_JC',         'HIP_R_IMU'}; ...    
%             {'body_rhip_rknee',     'KNEE_R_JC',        'KNEE_R_IMU'}; ...
%             {'body_rknee_rankle',   'ANKLE_R_JC',       'ANKLE_R_IMU'}...
            };
        
        % joint limits [deg]
        jointRangeOfMotion = {};
        jointRangeOfMotion_notes = {};
        
        % subject demographics (id, height [m], weight [kg])
        subjectDemographicArray = {...
            {1,  1.70, 65.7, 'm'} ...
            {2,  1.55, 70.0, 'f'} ...
            {3,  1.60, 47.7, 'f'} ...
            {4,  1.83, 90.0, 'm'} ...
            {5,  1.50, 52.0, 'f'} ...
            {6,  1.70, 70.0, 'm'} ...
            {7,  1.68, 71.0, 'm'} ...
            {8,  1.75, 60.0, 'm'} ...
            {9,  1.58, 52.0, 'f'} ...
            {10, 1.77, 63.5, 'm'} ...
            {11, 1.70, 61.0, 'm'} ...
            {12, 1.63, 60.8, 'f'} ...
            {13, 1.80, 70.0, 'm'} ...
            {14, 1.76, 68.0, 'm'} ...
            {15, 1.64, 63.0, 'f'} ...
            {16, 1.57, 55.0, 'f'} ...
            {17, 1.70, 82.0, 'm'} ...
            {18, 1.68, 56.0, 'f'} ...
            {19, 1.83, 71.0, 'm'} ...
            {20, 1.77, 75.0, 'm'} ...
            };
        
		linkLengthSymmetry = {...
			};
        trcSensorDecorator = {'position'};
        imuSensorDecorator = {'gyroscope', 'accelerometer'};
        imuSensorYaw = 0;
        
        lengthUseInd = 0;
    end
    
    methods
        function obj = rlModelInstance_healthy1(subjectId)
            obj = obj@rlModelInstance();
            
            obj.datasetName = 'healthy1';
            obj.subjectId = subjectId;
            obj.loadDemographics();
        end
        

        function  [dt, time, data] = loadData_trc(obj, pathToTrc, algorithmParam)
            pathToTrcHeader = strrep(pathToTrc, '.data', '.header');
            trc = loadDataFromTrc_healthy1_arsLoader(pathToTrcHeader, pathToTrc);
            
            dt = trc.dt;
            time = trc.time;
            data = trc.data;
        end
        
        function [data] = loadData_imu(obj, currFileEntry, algorithmParam)
            data = [ ...
                arsLoader(currFileEntry.filePathImuHip);...
                arsLoader(currFileEntry.filePathImuKneeRight);...
                arsLoader(currFileEntry.filePathImuAnkleRight);...
                ];
            sensorNames = {'Hip', 'KneeRight', 'KneeLeft'};

%            data = [ ...
%                 arsLoader(currFileEntry.filePathImuKneeRight);...
%                 arsLoader(currFileEntry.filePathImuAnkleRight);...
%                 ];
%             sensorNames = {'KneeRight', 'KneeLeft'};

            for i = 1:length(data)
                data(i).sensorId = data(i).name;
                data(i).name = sensorNames{i};
            end
        end

        
        function [kinematicTransform, dynamicTransform] = applyLinkTransform(obj, targetFrameStr, sourceFrameStr, linkVectorMean, lengthUseInd, dataInstance)
            linkLength = norm(linkVectorMean);
            
            % set any link-specific information
            switch targetFrameStr
                case {'length_rshoulder_rhip'} % either doesn't exist or already accounted for in other limbs
                    dumasFrameStr = 'torso';
                    
                case {'length_rhip_rknee'}
                    dumasFrameStr = 'thigh';
                    
                case {'length_rknee_rankle'}
                    dumasFrameStr = 'leg';
            end

            switch targetFrameStr
                case {'length_rshoulder_rhip', 'length_rhip_rknee', 'length_rknee_rankle'}
                    X00Length = linkVectorMean;
                    rotMtxDumas = rotz(pi/2)*rotx(pi/2);
            end
            
            switch obj.linkDefinition
                case 'X00'
                    correctedLinkLength = X00Length;
            end
            
            T = eye(4);
            T(1:3, 4) = correctedLinkLength;
            
            if ~isempty(dumasFrameStr)
                mass =            lookupTableDumas('mass',     dumasFrameStr, obj.gender, [], [])*obj.body_weight;
                comScale =        lookupTableDumas('com',      dumasFrameStr, obj.gender, linkLength, []);
                inertialScale =   lookupTableDumas('inertial', dumasFrameStr, obj.gender, linkLength, mass);
                
                % parallel axis, since Dumas inertia assumes it's on the edge of
                % the limb
%                 inertial_Hyugens =  mass * ...
%                     [comScale(2)^2+comScale(3)^2,  -comScale(1)*comScale(2),     -comScale(1)*comScale(3);
%                     -comScale(1)*comScale(2),       comScale(1)^2+comScale(3)^2, -comScale(2)*comScale(3);
%                     -comScale(1)*comScale(3),      -comScale(2)*comScale(3),      comScale(1)^2+comScale(2)^2];
                
                com = rotMtxDumas*comScale;
                inertial = rotMtxDumas*inertialScale;
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
                case {'ANKLE_R_MED' 'KNEE_R_MED'}
                    X00Length = [linkDistance 0 0]';
                    
                case {'ANKLE_R_LAT' 'KNEE_R_LAT'}
                    X00Length = [-linkDistance 0 0 ]';
                    
%                 case {'SHOULDER_L'}
%                     X00Length = [0 linkDistance 0]'; % something isn't right with the model...
                    
                case {'ASIS_L', 'SHOULDER_L'}
                    X00Length = [linkDistance 0 0]'; % something isn't right with the model...
                    
                case {'ASIS_R', 'SHOULDER_R'}
                    X00Length = [0 0 0]';
                    
                otherwise % IMU markers
                    X00Length = markerOffset;
            end
            
            switch obj.linkDefinition
                case 'X00'
                    correctedLinkLength = X00Length;
            end
            
            T = eye(4);
            T(1:3,4) = correctedLinkLength;
            
            sensorTransform = obj.assembleSenTransform(targetMarkerStr, sourceFrameStr, obj.trcSensorDecorator, T);
        end
        
        function imuTransforms = applyImuSensor(obj, R_sensor, imus, algorithmParam, dataInstance)
            switch algorithmParam.modalityCalibration
                case 'mocap'
                    imuTransforms = obj.calcImuAttachmentMocap(R_sensor, dataInstance);
            end
        end
        
        function obj = calcInit(obj)
            optInitParam.optTransXYZ = ones(3, 1);
            optInitParam.optRotXYZ = [0 0 1]';
            optInitParam.optJoints = zeros(length(obj.model.joints), 1);
            
            for i = 1:length(obj.model.joints)
                if strcmpi(obj.model.joints(i).name(end), '0') % sagittal joint
                    optInitParam.optJoints(i) = 1;
                end
            end
            
            optInitParam.markersIndTransXYZ = obj.findMarkerInd({'ANKLE_R_MED', 'ANKLE_R_LAT'});
            optInitParam.markersIndRotXYZ = obj.findMarkerInd({'KNEE_R_MED', 'KNEE_R_LAT'});
            
            [obj, T, pos] = calcOptInit(obj, optInitParam);
        end
    end
end