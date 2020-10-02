classdef rlModelInstance < handle % < rlModel
    properties
        % constants
        sensorNormThreshold = 100;
        
        % subject specific data
        datasetName = [];
        subjectId = [];
        body_height = [];
        body_weight = [];
        gender = '';
        
        % linkage attachment models
        linkDefinition = []; % 'InitialPose' 'X00'
        
%         sensorNameFull = [];
%         sensorAddStatus = [];
        
        % RL model
        model = [];
        model_ekf = [];
        
        % model templates for saving/loading
        kinematicTransform;
        dynamicTransform;
        sensorTransform;
        sensorSecondaryTransform;
        
        flatSensorAttachmentArray = {};
        flatSensorSecondaryAttachmentArray = {};
        
        % inds for tracking and errors
%        inds_mocapMarkers_rightArm;
%        inds_mocapMarkers_leftArm;
%        inds_mocapMarkers_rightLeg;
%        inds_mocapMarkers_leftLeg;
%        inds_mocapMarkers_torso;
%        inds_mocapMarkers;
        
%        inds_joints_rightArm;
%        inds_joints_leftArm;
%        inds_joints_rightLeg;
%        inds_joints_leftLeg;
%        inds_joints_torso;
%        inds_joints;
    end
    
    properties(Abstract = true)
         % link lengths (currentFrame, previousFrame, previousMarker, currentMarker)
        linkLengthArray;
        
        % marker arrangements (attachmentFrame, [markers])
        sensorAttachmentArray;
        
        % subject demographics (id, height [m], weight [kg])
        subjectDemographicArray;
    end
    
    methods(Abstract = true)
        obj = applyLinkTransform(obj, currLinkLengthArrangement, dataInstance);
        obj = applyMarkerSensor(obj, currLinkLengthArrangement, dataInstance);
    end
    
    methods
        function obj = rlModelInstance()
            obj.flatSensorAttachmentArray = rlModelInstance.flattenSensorAttachmentArray(obj.sensorAttachmentArray);
            obj.flatSensorSecondaryAttachmentArray = rlModelInstance.flattenSensorAttachmentArray(obj.sensorSecondaryAttachmentArray);
        end
        
        function initializeJointAndSensorInds(obj, strSensors, strJoints)
            % construct the sensor array if does not exist
            if ~exist('strSensor', 'var')
                indSensors = 0;
                strSensors = {};
                for i = 1:size(obj.sensorAttachmentArray, 1)
                    currSensorAttachmentArray = obj.sensorAttachmentArray{i};
                    for j = 3:size(currSensorAttachmentArray, 2)
                        indSensors = indSensors + 1;
                        strSensors{end+1} = currSensorAttachmentArray{j};
                    end
                end
            end
            
            % construct the joint array (assumes that a model has been loaded)
            if ~exist('strJoints', 'var') && ~isempty(obj.model)
                strJoints = {obj.model.joints.name};
            end
            
            % match the sensor strings
%             obj.inds_mocapMarkers_rightArm = find(ismember(strSensors, obj.mocapMarkers_rightArm)); 
%             obj.inds_mocapMarkers_leftArm = find(ismember(strSensors, obj.mocapMarkers_leftArm)); 
%             obj.inds_mocapMarkers_rightLeg = find(ismember(strSensors, obj.mocapMarkers_rightLeg)); 
%             obj.inds_mocapMarkers_leftLeg = find(ismember(strSensors, obj.mocapMarkers_leftLeg)); 
%             obj.inds_mocapMarkers_torso = find(ismember(strSensors, obj.mocapMarkers_torso)); 
%             obj.inds_mocapMarkers = sort([obj.inds_mocapMarkers_rightArm ...
%                 obj.inds_mocapMarkers_leftArm ...
%                 obj.inds_mocapMarkers_rightLeg ...
%                 obj.inds_mocapMarkers_leftLeg ...
%                 obj.inds_mocapMarkers_torso]);
%             
%             % match the sensor strings
%             obj.inds_joints_rightArm = find(ismember(strJoints, obj.joints_rightArm)); 
%             obj.inds_joints_leftArm = find(ismember(strJoints, obj.joints_leftArm)); 
%             obj.inds_joints_rightLeg = find(ismember(strJoints, obj.joints_rightLeg)); 
%             obj.inds_joints_leftLeg = find(ismember(strJoints, obj.joints_leftLeg)); 
%             obj.inds_joints_torso = find(ismember(strJoints, obj.joints_torso)); 
%             obj.inds_joints = sort([obj.inds_joints_rightArm ...
%                 obj.inds_joints_leftArm ...
%                 obj.inds_joints_rightLeg ...
%                 obj.inds_joints_leftLeg ...
%                 obj.inds_joints_torso]);
        end
        
        function [measurementOutData, indx] = rearrangeMeasurement(obj, measurementSourceData, measurementSourceLabel, targetLabelOrder, backupMeasurementInput)
            % rearrange input measurement source into the same order/len as
            % obj.sensorNameFull 
            
            if ~exist('targetLabelOrder', 'var') || isempty(targetLabelOrder)
                targetLabelOrder = obj.getSensorNameFull();
            end
            
            lenTime = size(measurementSourceData, 1);
            lenDof = length(targetLabelOrder);
            
            measurementOutData = SensorMeasurement(lenTime, lenDof);
            
            for i = 1:lenDof
                foundInd = find(ismember(measurementSourceLabel, targetLabelOrder{i}), 1);
                if ~isempty(foundInd)
                    indx(i) = foundInd;
                    measurementOutData(:, i) = measurementSourceData(:, indx(i));
                else
                    indx(i) = 0; % leave it empty (pull from the original, which will be full sized)
                    temp = backupMeasurementInput(:, i).getMesArray;
                    
                    measurementOutData(:, i) = backupMeasurementInput(:, i);  
                    measurementOutData(:, i).setMesArray(zeros(size(temp))); % but replace the entries with 0 to avoid confusion
                end
                
            end
        end
        
        function loadModel(obj, filepathModel)
            obj.model = rlCModel(filepathModel);
        end
        
        function sensorList = genSensorList(obj)
            indSensors = 0;
            for i = 1:size(obj.sensorAttachmentArray, 1)
                currSensorAttachmentArray = obj.sensorAttachmentArray{i};
                for j = 3:size(currSensorAttachmentArray, 2)
                    indSensors = indSensors + 1;
                    sensorList{indSensors} = currSensorAttachmentArray{j};
                end
            end
            
            for i = 1:size(obj.sensorSecondaryAttachmentArray, 1)
                currSensorAttachmentArray = obj.sensorSecondaryAttachmentArray{i};
                for j = 3:size(currSensorAttachmentArray, 2)
                    indSensors = indSensors + 1;
                    sensorList{indSensors} = currSensorAttachmentArray{j};
                end
            end
        end
        
        function [kinematicTransform, dynamicTransform, sensorTransform] = makeModel(obj, filepathModel, dataInstance, algorithmParam, linkLengthUseCell, sensorLengthUseCell)
            obj.model = rlCModel(filepathModel);
            obj.linkDefinition = algorithmParam.linkDefinition;
            
            if ~exist('linkLengthUseCell', 'var')
                lengthUseIndSingle = obj.lengthUseInd;
                if lengthUseIndSingle == 0
                    lengthUseIndSingle = 1:length(dataInstance.time);
                end
                
                for i = 1:size(obj.linkLengthArray, 1)
                    linkLengthUseCell{i} = lengthUseIndSingle;
                end
                for i = 1:size(obj.flatSensorAttachmentArray, 1)
                    sensorLengthUseCell{i} = lengthUseIndSingle;
                end
            end
           
            % calculate link offsets
            for i = 1:size(obj.linkLengthArray, 1)
                % load link length
                [sourceFrameStr, targetFrameStr, linkLength, linkVectorMean, linkVectorStd] = ...
                    rlModelInstance.linkLengthArrayData(obj.linkLengthArray{i}, linkLengthUseCell{i}, dataInstance);
                
                % figure out how to attach it
                [kinematicTransform(i), dynamicTransform(i)] = ...
                    obj.applyLinkTransform(targetFrameStr, sourceFrameStr, linkVectorMean, linkLengthUseCell{i}, dataInstance);
                
                fprintf('rlModelInstance: calculating length %s: [%f, %f, %f] ... ', targetFrameStr, ...
                    kinematicTransform(i).t(1, 4), kinematicTransform(i).t(2, 4), kinematicTransform(i).t(3, 4));
                fprintf('applying %f \n ', norm(kinematicTransform(i).t(1:3, 4)));
                
                % if the link length is malformed, NaN it for later correction
                if isnan(linkLength)
                    kinematicTransform(i).t(1:3, 4) = NaN; % blank out the joint lengths
                end
            end 
            
            % if there is any NAN entries, check the right and left to see if
            % the link lengths can be mirrored
            [kinematicTransform, dynamicTransform] = obj.fillInSym('kinematic', ...
                kinematicTransform, dynamicTransform, obj.linkLengthSymmetry, linkLengthUseCell, dataInstance);

            % now add as transform
            obj.addKinDynTransformToModel(obj.model, kinematicTransform, dynamicTransform);
            obj.model.forwardPosition();
            
            % add sensors
            indSensors = size(obj.flatSensorAttachmentArray, 1);
%             obj.sensorAddStatus = [];
%             obj.sensorNameFull = {};
            for i = 1:size(obj.flatSensorAttachmentArray, 1)
                % load marker
                [sourceFrameStr, targetMarkerStr, linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData] = ...
                    rlModelInstance.sensorLengthArrayData(obj.flatSensorAttachmentArray{i}, sensorLengthUseCell{i}, dataInstance);
                
                % attach marker
                checkSensBool = checkSensorPlacement(obj.model, targetMarkerData, linkVectorMean, ...
                    sensorLengthUseCell{i}, targetMarkerStr, sourceFrameStr, obj.sensorNormThreshold);
                
                sensorTransformTmp = obj.applyMarkerSensor(sourceFrameStr, targetMarkerStr, ...
                    linkVectorMean, sensorLengthUseCell{i}, dataInstance);
                
                if checkSensBool
                    fprintf('rlModelInstance: attaching sensor (%u) %s to frame %s: %f\n', ...
                        i, targetMarkerStr, sourceFrameStr, linkLength);
                    sensorTransform(i) = sensorTransformTmp;
                    %                         obj.sensorAddStatus(indSensors) = 1;
                else
                    fprintf('rlModelInstance: blanking sensor (%u) %s to frame %s: %f\n', ...
                        i, targetMarkerStr, sourceFrameStr, linkLength);
                    sensorTransform(i) = obj.assembleSensorTransformEmpty(targetMarkerStr, sourceFrameStr, sensorTransformTmp.decorators);
                    %                         obj.sensorAddStatus(indSensors) = 0;
                end
            end
            
            % if there is any NAN entries, check the right and left to see if
            % the link lengths can be mirrored
%            [sensorTransform, ~] = obj.fillInSym('sensor', ...
%                 sensorTransform, [], obj.sensorSymmetry, lengthUseInd, dataInstance, currSensorAttachment);
            
            obj.addSenTransformToModel(obj.model, sensorTransform);
            
            obj.kinematicTransform = kinematicTransform;
            obj.dynamicTransform = dynamicTransform;
            obj.sensorTransform = sensorTransform;
            
%             if length(obj.sensorAddStatus) < length(dataInstance.measurement_labels)
%                 % there's more entries in measurement_label, probably because
%                 % unknownmarkers are added via
%                 % algorithmParam.addUnknownMarkersToEkf
%                 
%                 for i = (length(obj.sensorAddStatus)+1):length(dataInstance.measurement_labels)
%                     % if there are issues with the data (ie 
%                     if dataInstance.measurement(1, i).type == 0
%                         obj.sensorAddStatus(i) = 0;
%                     else
%                         obj.sensorAddStatus(i) = 1;
%                     end
%                 end
%             end
            
%             % remove missing markers
%             if min(obj.sensorAddStatus) == 0
%                 % the marker was not attached. drop it from the data mes
%                 measurement_tmp = dataInstance.measurement(:, obj.sensorAddStatus == 1);
%                 measurementLabel_tmp = dataInstance.measurement_labels(:, obj.sensorAddStatus == 1);
%                 measurementSensor_tmp = dataInstance.measurement_sensor(:, obj.sensorAddStatus == 1);
%                 
%                 dataInstance.measurement = measurement_tmp;
%                 dataInstance.measurement_labels = measurementLabel_tmp;
%                 dataInstance.measurement_sensor = measurementSensor_tmp;
%             end
        end
        
        function runningLenVector = symVal(obj, runMode, mainTransform, secondaryTransform, symmetryList, i)
            symInd = [];
            % find the corresponding symmetry array
            for j = 1:length(symmetryList)
                switch runMode
                    case 'kinematic'
                        symInd = find(ismember(symmetryList{j}, mainTransform(i).frameName));
                        
                    case 'sensor'
                        symInd = find(ismember(symmetryList{j}, mainTransform(i).sensorName));
                end
                
                if ~isempty(symInd)
                    break;
                end
            end
            
            % once we've found the proper value, calculate the average link
            % length according to non-NAN entries in the symmetry list
            if ~isempty(symInd)
                switch runMode
                    case 'kinematic'
                        transformNames = {mainTransform.frameName};
                        
                    case 'sensor'
                        transformNames = {mainTransform.sensorName};
                end
                
                % create a non-nan average length
                runningLenVector = zeros(3, 1);
                divFactor = 0;
                symInd = find(ismember(transformNames, symmetryList{j}));
                for k = 1:length(symmetryList{j})
                    if isempty(mainTransform(symInd(k)).t)
                        continue
                    end
                    symLenVector = mainTransform(symInd(k)).t(1:3, 4);
                    
                    if isempty(find(isnan(symLenVector)))
                        divFactor = divFactor + 1;
                        runningLenVector = runningLenVector + abs(symLenVector);
                    end
                end
                
                if divFactor > 0
                    runningLenVector = runningLenVector/divFactor;
                else 
                    runningLenVector = [];
                end
            else
                runningLenVector = [];
            end
        end

        function [mainTransformOut, secondaryTransformOut] = fillInSym(obj, runMode, ...
                mainTransform, secondaryTransform, symmetryList, lengthUseInd, dataInstance, currSensorAttachment)

            % make a copy for exporting
            mainTransformOut = mainTransform;
            secondaryTransformOut = secondaryTransform;
            
            for i = 1:length(mainTransform)
                if ~isempty(find(isnan(mainTransform(i).t), 1)) || isempty(mainTransform(i).t)
                    % find if this entry has a corresponding on the other side
                    runningLenVector = obj.symVal(runMode, mainTransform, secondaryTransform, symmetryList, i);
                    
                    if ~isempty(runningLenVector)
                        switch runMode
                            case 'kinematic'
                                fprintf('rlModelInstance: re-calculating length %s ... ', mainTransform(i).frameName);
                                
                                [sourceFrameStr, targetFrameStr, linkLength, linkVectorMean, linkVectorStd] = ...
                                    rlModelInstance.linkLengthArrayData(obj.linkLengthArray{i}, lengthUseInd{i}, dataInstance);
                                linkLength = runningLenVector;
                                linkDistance = norm(linkLength);
                                
                                % targetFrameStr, sourceFrameStr, linkVectorMean, linkLengthUseCell{i}, dataInstance
                                [mainTransformOut(i), secondaryTransformOut(i)] = obj.applyLinkTransform(targetFrameStr, sourceFrameStr, linkLength, lengthUseInd{i}, dataInstance);
                                
                            case 'sensor'
                                fprintf('rlModelInstance: re-calculating length %s ... ', mainTransform(i).sensorName);
                                
                                [sourceFrameStr, targetMarkerStr, linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData] = ...
                                    obj.sensorLengthArrayData(obj.flatSensorAttachmentArray{i}, sensorLengthUseCell{i}, dataInstance);
                                linkLength = runningLenVector;
                                linkDistance = norm(linkLength);
                                
                                % attach marker
                                checkSensBool = checkSensorPlacement(obj.model, targetMarkerData, linkVectorMean, ...
                                    sensorLengthUseCell{i}, targetMarkerStr, sourceFrameStr, obj.sensorNormThreshold);
                                
                                sensorTransformTmp = obj.applyMarkerSensor(sourceFrameStr, targetMarkerStr, ...
                                    linkVectorMean, sensorLengthUseCell{i}, dataInstance);
                                
                                if checkSensBool
                                    fprintf('rlModelInstance: attaching sensor (%u) %s to frame %s: %f\n', ...
                                        i, targetMarkerStr, sourceFrameStr, linkLength);
                                    mainTransformOut(i) = sensorTransformTmp;
                                    %                         obj.sensorAddStatus(indSensors) = 1;
                                else
                                    fprintf('rlModelInstance: blanking sensor (%u) %s to frame %s: %f\n', ...
                                        i, targetMarkerStr, sourceFrameStr, linkLength);
                                    mainTransformOut(i) = obj.assembleSensorTransformEmpty(targetMarkerStr, sourceFrameStr, sensorTransformTmp.decorators);
                                    %                         obj.sensorAddStatus(indSensors) = 0;
                                end
                        end
                        
                        fprintf('applying %f \n ', norm(mainTransformOut(i).t(1:3, 4)));
                    end
                        
                    % however, if there is still a problem with the
                    % transform after this match, we will need to reset it
                    if ~isempty(find(isnan(mainTransform(i).t), 1)) || isempty(mainTransform(i).t)
                        switch runMode
                            case 'kinematic'
                                fprintf('rlModelInstance: better link values for %s not found from symmetry \n ', mainTransform(i).frameName);
                                mainTransformOut(i).t = eye(4);
                                
                                secondaryTransformOut(i).m = 0;
                                secondaryTransformOut(i).com = zeros(3, 1);
                                secondaryTransformOut(i).I = eye(3);
                                
                            case 'sensor'
                                fprintf('rlModelInstance: better link values for %s not found from symmetry \n ', mainTransform(i).sensorName);
                        end
                    end
                else
                    % keep the orig intact
                end
            end
        end
        
        function plotCurrentPose(obj, mdl)
            % plot the array of used markers
            %                 offsetArray = [dataInstance.data.HIP_BASE zeros(size(dataInstance.data.HIP_BASE))];
            %                 for i = 1:length(obj.measurement_indArray)
            %                     obj.measurement_indArray{i} = obj.measurement_indArray{i} - offsetArray;
            %                 end
            %                 for i = 1:length(obj.model.bodies)
            %                     vis.addMarker(['COM_' obj.model.bodies(i).name], ...
            %                         obj.model.bodies(i).t(1:3, 4) + obj.model.bodies(i).com, [0.8 0 0 1]);
            %                 end
            %                 testModelAllLinks(obj.model);
            
            % plot the full set of raw data
            
            if nargin == 1
                if ~isempty(obj.model_ekf)
                    mdl = obj.model_ekf;
                else
                    mdl = obj.model;
                end
            end
            
            vis = rlVisualizer('vis',640,480); 
            mdl.forwardPosition();
            vis.addModel(mdl);
            
            applyMarkersToVisualization(vis, mdl, [], [], [], [])
            
            vis.update();
%             pause();
           
%                 testModelAllLinks(obj.model);
        end   
        
        function modelStruct = saveModel(obj, filepathModelParam)
            modelStruct.kinematicTransform = obj.kinematicTransform;
            modelStruct.dynamicTransform = obj.dynamicTransform;
            modelStruct.sensorTransform = obj.sensorTransform;
            modelStruct.sensorSecondaryTransform = obj.sensorSecondaryTransform;
            
            if ~isempty(filepathModelParam)
                save(filepathModelParam, 'modelStruct');
            end
        end

        function loadTrcEkfModelFromModelSpecs(obj, filepathModel, filepathModelParam, sensorList)
            load(filepathModelParam);    
            
            obj.model_ekf = rlCModel(filepathModel);
            if isempty(sensorList)
                sensorList = obj.genSensorList();
            end
            
            obj.loadModelFromModelStruct(obj.model_ekf, modelStruct, sensorList);
            
            if ~isempty(modelStruct.sensorSecondaryTransform)
                obj.addSecondaryTransformToModel(obj.model_ekf, modelStruct.sensorSecondaryTransform, sensorList);
            end
            obj.sensorSecondaryTransform = modelStruct.sensorSecondaryTransform;
        end
        
        function loadImuEkfModelFromModelSpecs(obj, filepathModel, filepathModelParam, secondarySensorList)
            load(filepathModelParam);
            
            obj.model_ekf = rlCModel(filepathModel);
            modelStruct.sensorTransform = []; % remove trc models
            sensorList = [];
            obj.loadModelFromModelStruct(obj.model_ekf, modelStruct, sensorList);

%             obj.model_ekf = obj.copyJointPosition(obj.model, obj.model_ekf);
            
            sensorSecondaryTransform = modelStruct.sensorSecondaryTransform;
            if ~isempty(sensorSecondaryTransform)
                for i = 1:length(sensorSecondaryTransform)
                    if strcmpi(sensorSecondaryTransform(i).sensorName(end-3:end), '_YAW')
                        % replace the loaded yaw sensor with a current one
                        frame = obj.model_ekf.getFrameByName(sensorSecondaryTransform(i).frameName);
                        T = eye(4);
                        T(1:3, 1:3) = frame.t(1:3,1:3)';
                        T(1:3, 4) = zeros(3, 1);
                        sensorSecondaryTransform(i).t = T;
                    end
                end

                if exist('secondarySensorList', 'var')
                    keepInds = obj.addSecondaryTransformToModel(obj.model_ekf, sensorSecondaryTransform, secondarySensorList);
                else
                    keepInds = obj.addSecondaryTransformToModel(obj.model_ekf, sensorSecondaryTransform);
                end
                
            end
            
            obj.sensorSecondaryTransform = sensorSecondaryTransform(keepInds);
        end
        
        function loadModelFromModelSpecs(obj, filepathModel, filepathModelParam, sensorList)
            load(filepathModelParam);
            
            if isempty(sensorList)
                sensorList = obj.genSensorList();
            end
            
            if exist('saveVar', 'var')
                modelStruct = saveVar.modelStruct;
            end
            
            obj.model = rlCModel(filepathModel);
            obj.loadModelFromModelStruct(obj.model, modelStruct, sensorList);
        end
        
        function saveVar = loadModelFromModelSpecsNoSensor(obj, filepathModel, filepathModelParam, updateDynamicModel)
            load(filepathModelParam);
            
            if ~exist('updateDynamicModel', 'var')
                updateDynamicModel = 0;
            end
            
            if updateDynamicModel 
                % there's no dynamic model properly passed in, will update
                obj.linkDefinition = 'initialpose';
                saveVar.modelStruct.dynamicTransform = generateDynamicTransform(obj, saveVar.modelStruct);
            end
            
            obj.model = rlCModel(filepathModel);
            sensorList = [];
            obj.loadModelFromModelStruct(obj.model, saveVar.modelStruct, sensorList);
        end
        
        function dynamicTransform = generateDynamicTransform(obj, modelStruct)
            % calculate link offsets
            for i = 1:numel(modelStruct.kinematicTransform)
                % load link length
                %                 [sourceFrameStr, targetFrameStr, linkLength, linkVectorMean, linkVectorStd] = ...
                %                     rlModelInstance.linkLengthArrayData(obj.linkLengthArray{i}, linkLengthUseCell{i}, dataInstance);
                targetFrameStr = modelStruct.kinematicTransform(i).frameName;
                sourceFrameStr = obj.linkLengthArray{i}{2};
                linkVectorMean = modelStruct.kinematicTransform(i).t(1:3, 4);
                linkLengthUseCell = [];
                dataInstance = [];
                
                % figure out how to attach it
                [~, dynamicTransform(i)] = ...
                    obj.applyLinkTransform(targetFrameStr, sourceFrameStr, linkVectorMean, linkLengthUseCell, dataInstance);
                
            end
        end
        
        function loadModelFromModelStruct(obj, mdl, modelStruct, sensorList)
            % define links
            obj.addKinDynTransformToModel(mdl, modelStruct.kinematicTransform, modelStruct.dynamicTransform);
            obj.addSenTransformToModel(mdl, modelStruct.sensorTransform, sensorList);
            
            obj.kinematicTransform = modelStruct.kinematicTransform;
            obj.dynamicTransform = modelStruct.dynamicTransform;
            obj.sensorTransform = modelStruct.sensorTransform;
        end
        
        function targetMdl = copyTransforms(obj, sourceMdl, targetMdl)
            transform_names = {sourceMdl.transforms.name};
            for i=1:numel(targetMdl.transforms)
                currTargetFrameName = targetMdl.transforms(i).name;
                indx = find(ismember(transform_names,currTargetFrameName),1);
                if ~isempty(indx)
                    targetMdl.transforms(i).t =...
                        sourceMdl.transforms(indx).t;
                else
                    lala = 1;
                end
            end
        end
        
        function targetMdl = copyJointPosition(obj, sourceMdl, targetMdl)
            joint_names = {sourceMdl.joints.name};
            for i=1:numel(targetMdl.joints)
                indx = find(ismember(joint_names,targetMdl.joints(i).name),1);
                if ~isempty(indx)
                    targetMdl.joints(i).position =...
                        sourceMdl.joints(indx).position;
                end
            end
            targetMdl.forwardPosition();
        end
        
        function targetMdl = copySensors(obj, sourceMdl, targetMdl, sensorNamesToCopy)
            if nargin == 3
                sensorNamesToCopy = {};
            end
            
            transform_names = {sourceMdl.transforms.name};
            for i=1:numel(sourceMdl.sensors)
                indxTrans = find(ismember(transform_names, [sourceMdl.sensors(i).name 't']),1);
                indxToCopy = find(ismember(sensorNamesToCopy, sourceMdl.sensors(i).name));
                
                if ~isempty(sensorNamesToCopy) && isempty(indxToCopy)
                    % there is a sensornametocopy array, but it doesn't hold the
                    % proper value
                    continue;
                end
                
                if ~isempty(indxTrans)
                    M1 = SensorCore(sourceMdl.sensors(i).name);
                    typeArray = strsplit(sourceMdl.sensors(i).type, ',');
                    for j = 1:length(typeArray)-1
                        M1.addDecorator(typeArray{j});
                    end
                    
                    T = sourceMdl.transforms(indxTrans).t;
                    frameIn = sourceMdl.transforms(indxTrans).frame_in.name;
                    targetMdl.addSensor(M1, frameIn, T);
                end
            end
        end 
        
        function sensorSecondaryTransform = point2RotMethod(obj, R_sensor, dataInstance, jointAngles, lengthUseCell)
            [transform, sensorLengthUseCell] = obj.calc3MarkerAttachmentMocap(R_sensor, dataInstance, jointAngles, lengthUseCell);
                
            indSensors = length(obj.model_ekf.sensors);
            indSecondarySensor = 0;
            for i = 1:length(transform)
                % check the quality of the reconstruction and add if appropriate
                if ~isempty(transform(i).T_f_imu)
                    fprintf('rlModelInstance: attaching sensors %s to frame %s\n', ...
                        obj.sensorSecondaryAttachmentArray{i}{3}, ...
                        transform(i).frame);
                    
                    currSensorAttachmentArray = obj.sensorSecondaryAttachmentArray{i};
                    for j = 3:length(currSensorAttachmentArray)
                        % determine new offset
                        currSensorAttachment = [currSensorAttachmentArray(1) currSensorAttachmentArray(2) currSensorAttachmentArray(j)];
                        
                        [sourceFrameStr, targetMarkerStr, linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData, useIntersect] = ...
                            rlModelInstance.sensorLengthArrayData(currSensorAttachment, sensorLengthUseCell{i}, dataInstance);

                        targetMarkerDataLen = mean(targetMarkerData(sensorLengthUseCell{i}, :), 1);
                        T = transform(i).T_f_imu;
                        T(1:3, 4) = targetMarkerDataLen';
                        T_out = transform(i).T_f_0*T;
                        
                        indSecondarySensor = indSecondarySensor + 1;
                        indSensors = indSensors + 1;
                        sensorSecondaryTransform(indSecondarySensor) = obj.assembleSenTransform3dTrc(...
                            obj.sensorSecondaryAttachmentArray{i}{j}, transform(i).frame, obj.trcSensorDecorator, T_out(1:3, 1:3), T_out(1:3, 4));
                        %                         obj.sensorAddStatus(indSensors) = 1;
                    end
                else
                    fprintf('rlModelInstance: blanking sensors %s to frame %s\n', ...
                        obj.sensorSecondaryAttachmentArray{i}{3}, ...
                        transform(i).frame);
                    
                     for j = 3:length(obj.sensorSecondaryAttachmentArray{i})
                         indSecondarySensor = indSecondarySensor + 1;
                         indSensors = indSensors + 1;
                         sensorSecondaryTransform(indSecondarySensor) = obj.assembleSensorTransformEmpty(...
                             obj.sensorSecondaryAttachmentArray{i}{j}, transform(i).frame, obj.trcSensorDecorator);               
                     end
                end
            end
        end
        
        function sensorSecondaryTransform = initialPoseMethod(obj, R_sensor, dataInstance, jointAngles, lengthUseCell)
            % add sensors
            for i = 1:size(obj.flatSensorSecondaryAttachmentArray, 1)
                sensorLengthUseCell{i} = 1:10;
                
                % load marker
                [sourceFrameStr, targetMarkerStr, linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData] = ...
                    rlModelInstance.sensorLengthArrayData(obj.flatSensorSecondaryAttachmentArray{i}, sensorLengthUseCell{i}, dataInstance);
                
                % attach marker
                checkSensBool = checkSensorPlacement(obj.model, targetMarkerData, linkVectorMean, ...
                    sensorLengthUseCell{i}, targetMarkerStr, sourceFrameStr, obj.sensorNormThreshold);
                
                sensorTransformTmp = obj.applyMarkerSensor(sourceFrameStr, targetMarkerStr, ...
                    linkVectorMean, sensorLengthUseCell{i}, dataInstance);
                
                if checkSensBool
                    fprintf('rlModelInstance: attaching sensor (%u) %s to frame %s: %f\n', ...
                        i, targetMarkerStr, sourceFrameStr, linkLength);
                    sensorSecondaryTransform(i) = sensorTransformTmp;
                    %                         obj.sensorAddStatus(indSensors) = 1;
                else
                    fprintf('rlModelInstance: blanking sensor (%u) %s to frame %s: %f\n', ...
                        i, targetMarkerStr, sourceFrameStr, linkLength);
                    sensorSecondaryTransform(i) = obj.assembleSensorTransformEmpty(targetMarkerStr, sourceFrameStr, sensorTransformTmp.decorators);
                    %                         obj.sensorAddStatus(indSensors) = 0;
                end
            end
        end
        
        function sensorSecondaryTransform = makeTRCModel(obj, filepathModel, dataInstance, initPose, jointAngles, algorithmParam, lengthUseCell)
            %Load xml Model
            jointAngles = obj.makeEkfModel(filepathModel, dataInstance, initPose, jointAngles, algorithmParam);
            
            % now copy over the base sensors
            obj.model_ekf = obj.copySensors(obj.model, obj.model_ekf);
            
            % figure out the transforms for the non-X00 sensors and attach them
            R_sensor = eye(3);
            if ~isempty(obj.sensorSecondaryAttachmentArray)
                switch obj.linkDefinition
                    case 'initialpose'
                        sensorSecondaryTransform = obj.initialPoseMethod(R_sensor, dataInstance, jointAngles, lengthUseCell);
                        
                    case 'X00'
                        sensorSecondaryTransform = obj.point2RotMethod(R_sensor, dataInstance, jointAngles, lengthUseCell);
                end
                
                obj.sensorSecondaryTransform = sensorSecondaryTransform;
                obj.addSecondaryTransformToModel(obj.model_ekf, obj.sensorSecondaryTransform);     
            else
                sensorSecondaryTransform = [];
            end
            
            obj.model_ekf.forwardPosition();
        end
        
        function jointAngles = makeEkfModel(obj, filepathModel, dataInstance, initPose, jointAngles, algorithmParam)
            obj.model_ekf = rlCModel(filepathModel);

            % set the model to be in the first pose
            obj.model.position = initPose;
            obj.model.base = algorithmParam.baseFrame;
            obj.model.forwardPosition();
            
            jointRemapping = contains({obj.model.joints.name}, {obj.model_ekf.joints.name});
%             initPose = initPose(:, jointRemapping);
            jointAngles = jointAngles(:, jointRemapping);
            
            obj.model_ekf.base = algorithmParam.baseFrame;
            
            %We copy all the transforms for the IMU model from the Mocap
            %based model
            obj.model_ekf = obj.copyTransforms(obj.model, obj.model_ekf);
            
            %Here we set initial position
            obj.model_ekf = obj.copyJointPosition(obj.model, obj.model_ekf);
        end
        
        function sensorSecondaryTransform = makeIMUModel(obj, filepathModel, dataInstance, R_sensor, initPose, jointAngles, algorithmParam, lengthUseCell)
            %Load xml Model
            jointAngles = obj.makeEkfModel(filepathModel, dataInstance, initPose, jointAngles, algorithmParam);
            
            %At this point we have identically located fixed and moving
            %base models, now attach IMU sensors to the IMU model 
            % add sensors
            transform = obj.applyImuSensor(R_sensor, algorithmParam, jointAngles, dataInstance, lengthUseCell);
     
            indSensors = length(obj.model_ekf.sensors);
            indSecondarySensor = 0;
            for i = 1:length(transform)
                % update dataInstance.measurement_sensor (flip 0 to 1)
                sensorMarkerStr = [transform(i).name];
%                 obj.sensorNameFull{end+1} = sensorMarkerStr;
                imuDecorator = obj.imuSensorDecorator;
                
                % check the quality of the reconstruction and add if appropriate
                if ~isempty(transform(i).R)
                    fprintf('rlModelInstance: attaching sensors %s to frame %s\n', sensorMarkerStr, transform(i).frame);
                    
                    % update dataInstance.sensorAddStatus (append 0 or 1)
                    indSecondarySensor = indSecondarySensor + 1;
                    indSensors = indSensors + 1;
                    sensorSecondaryTransform(indSecondarySensor) = obj.assembleSenTransform3dTrc(...
                        sensorMarkerStr, transform(i).frame, imuDecorator, transform(i).R, transform(i).t);
%                     obj.sensorAddStatus(indSensors) = 1;
                else
                    fprintf('rlModelInstance: blanking sensors %s to frame %s\n', sensorMarkerStr, transform(i).frame);
                    
                    indSecondarySensor = indSecondarySensor + 1;
                    indSensors = indSensors + 1;
                    sensorSecondaryTransform(indSecondarySensor) = obj.assembleSensorTransformEmpty(...
                        sensorMarkerStr, transform(i).frame, imuDecorator);
%                     obj.sensorAddStatus(indSensors) = 0;
                end
            end
            
            if obj.imuSensorYaw
                imuDecorator = {};
                
                if isprop(obj, 'sensorSecondaryImuAttachmentArray')
                    localSensorAttachmentArray = obj.sensorSecondaryImuAttachmentArray;
                else
                    localSensorAttachmentArray = obj.sensorSecondaryAttachmentArray;
                end
                
                if isprop(obj, 'sensorYawSecondaryAttachmentArray')
                    sensorYawUse = obj.sensorYawSecondaryAttachmentArray;
                else
                    sensorYawUse = obj.sensorSecondaryAttachmentArray;
                end
                
                for i = 1:length(sensorYawUse)
                    ind = [];
                    currYawUseName = sensorYawUse{i}{2};
                    for j = 1:length(transform)
                        % find the matching transform
                        if strcmpi(localSensorAttachmentArray{j}{2}, currYawUseName)
                            ind = j;
                            break;
                        end
                    end
                    
                    sensorMarkerStr = [transform(ind).name '_YAW'];
%                     obj.sensorNameFull{end+1} = sensorMarkerStr;
                    imuDecorator{1} = 'yaw';
                    
                    frame = obj.model_ekf.getFrameByName(transform(ind).frame);
                    R = frame.t(1:3,1:3)';
                    t = zeros(3, 1);
                    
                    fprintf('rlModelInstance: attaching sensors %s to frame %s\n', ...
                        sensorMarkerStr, transform(ind).frame);
                    
                    % update dataInstance.sensorAddStatus (append 0 or 1)
                    indSecondarySensor = indSecondarySensor + 1;
                    indSensors = indSensors + 1;
                    sensorSecondaryTransform(indSecondarySensor) = obj.assembleSenTransform3dTrc(...
                        sensorMarkerStr, transform(ind).frame, imuDecorator, R, t);
%                     obj.sensorAddStatus(indSensors) = 1;
                end
            end
            
            obj.sensorSecondaryTransform = sensorSecondaryTransform;
            obj.addSecondaryTransformToModel(obj.model_ekf, obj.sensorSecondaryTransform);
            
%             for i = 1:length(transform)
%                 sensorMarkerStr = [transform(i).name ''];
%                 attachmentFrameStr = transform(i).frame;
%                 
%                 sensorNameFull{end+1} = sensorMarkerStr;
%                 sensorAddStatus(end+1) = 1;
%                 
%                 T = eye(4);
%                 T(1:3, 1:3) = transform(i).R;
%                 T(1:3, 4) = transform(i).t;
%                 
%                 M1 = SensorCore(sensorMarkerStr);
%                 for j = 1:length(obj.imuSensorDecorator)
%                     M1.addDecorator(obj.imuSensorDecorator{j});
%                 end
%                 obj.model_ekf.addSensor(M1, attachmentFrameStr, T);
%             end
%             
%             %Should we add yaw sensors or not?
%             if obj.imuSensorYaw
%                 for i = 1:size(obj.sensorImuAttachmentArray, 1)
%                     sensorMarkerStr = [transform(i).name '_YAW'];
%                     attachmentFrameStr = transform(i).frame;
%                     
%                     sensorNameFull{end+1} = sensorMarkerStr;
%                     sensorAddStatus(end+1) = 1;
%                     
%                     frame = obj.model_ekf.getFrameByName(attachmentFrameStr);
%                     T = eye(4);
%                     T(1:3,1:3) = frame.t(1:3,1:3)';
%                     
%                     sens = SensorCore(sensorMarkerStr);
%                     sens.addDecorator('yaw');
%                     obj.model_ekf.addSensor(sens,attachmentFrameStr,T);
%                 end
%             end
%             
%             % update the sensorNameFull and sensorAddStatus
%             obj.sensorNameFull = sensorNameFull;
%             obj.sensorAddStatus = sensorAddStatus;
%             obj.model_ekf.forwardPosition();

            if 0
                indArray = 1;
                vis = rlVisualizer('vis',640,480); obj.model_ekf.forwardPosition(); vis.addModel(obj.model_ekf); vis.update();
                applyMarkersToVisualization(vis, obj.model_ekf, dataInstance, indArray, [], []);
            end
        end
        
        function [transform, sensorLengthUseCell] = calc3MarkerAttachmentMocap(obj, R_sensor, dataInstance, jointAngles, sensorLengthUseCell)
            if ~exist('sensorLengthUseCell', 'var') || isempty(sensorLengthUseCell)
                lengthUseIndSingle = obj.lengthUseInd;
                if lengthUseIndSingle == 0
					% lengthUseIndSingle = 1:length(dataInstance.time);                
                    lengthUseIndSingle = 1;
                end
                
                for i = 1:size(obj.sensorSecondaryAttachmentArray, 1)
                    sensorLengthUseCell{i} = lengthUseIndSingle;
                end
            end
            
%             for i = 1:length(obj.flatSensorSecondaryAttachmentArray)
%                 flatSensorSecondaryLastCell(i) = obj.flatSensorSecondaryAttachmentArray{i}(end);
%             end

            if isprop(obj, 'sensorSecondaryImuAttachmentArray')
                localSensorAttachmentArray = obj.sensorSecondaryImuAttachmentArray;
            else
                localSensorAttachmentArray = obj.sensorSecondaryAttachmentArray;
            end
            
            addMarker = zeros(size(localSensorAttachmentArray, 1), 3);
            expectedMarkers = zeros(size(localSensorAttachmentArray, 1), 1);
            for i = 1:size(localSensorAttachmentArray, 1)
                currSensorAttachmentArray = localSensorAttachmentArray{i};
                expectedMarkers(i) = size(currSensorAttachmentArray, 2) - 2;
                
                for j = 3:size(currSensorAttachmentArray, 2)
                    currSensorAttachment = [currSensorAttachmentArray(1) currSensorAttachmentArray(2) currSensorAttachmentArray(j)];
                    
%                     indFlatSensor = find(ismember(flatSensorSecondaryLastCell, currSensorAttachment{3}));
                    indFlatSensor = i;
                    
                    [sourceFrameStr, targetMarkerStr, linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData, useIntersect] = ...
                        rlModelInstance.sensorLengthArrayData(currSensorAttachment, sensorLengthUseCell{indFlatSensor}, dataInstance);
                    
                    % attach marker
                    addMarker(i, j-2) = checkSensorPlacement(obj.model_ekf, targetMarkerData, linkVectorMean, ...
                        sensorLengthUseCell{indFlatSensor}, targetMarkerStr, sourceFrameStr, obj.sensorNormThreshold);

                    switch j
                        case 3
                            P_ind{i} = useIntersect;
                            P_data{i} = targetMarkerData;
                            
                        case 4
                            Q_ind{i} = useIntersect;
                            Q_data{i} = targetMarkerData;
                            
                        case 5
                            R_ind{i} = useIntersect;
                            R_data{i} = targetMarkerData;
                    end
                end
            end
            
            if ~iscell(R_sensor)
                R_sensor_new = repmat({R_sensor}, size(localSensorAttachmentArray, 1), 1);
                R_sensor = R_sensor_new;
            end
            
            for i = 1:size(localSensorAttachmentArray, 1)
                currSensorAttachmentArray = localSensorAttachmentArray{i};
                %Get the frame we are attaching the sensor to
                frame = obj.model_ekf.getFrameByName(currSensorAttachmentArray{1});
                
                T_f_imu_raw = {};
                T_f_0_avg = {};
                R_f_imu = {};
                
                for j = 1:size(jointAngles, 1)
                    T_f_imu_raw{j} = [];
                    T_f_0_avg{j} = [];
                    R_f_imu{j} = [];
                    
                    % set the model pose to the IK
                    obj.model_ekf.position = jointAngles(j, :);
                    obj.model_ekf.forwardPosition();
                    
                    %Transformation from frame to world
                    T_0_f = frame.t;
                    %Transformation from world to frame
                    T_f_0 = T_0_f;
                    T_f_0(1:3,1:3) = T_f_0(1:3,1:3)';
                    T_f_0(1:3,4) = -T_f_0(1:3,1:3)*T_f_0(1:3,4);
                    %Figure out transformation from sensor to world if the
                    %there are 3 markers on the sensor. The axes of the sensor
                    %will be defined as follows:
                    %X: Marker 1 -> Marker 2
                    %Z: Cross(X, Marker 1 -> Marker 3)
                    %Y: Cross(Z,X)
%                     R_0_imu{j} = T_f_0(1:3,1:3)'*R_sensor{i};       %Rotation Matrix IMU -> World
%                     %Position of IMU in global frame
%                     t_0_imu{j} = [];

                    clearCurrInd = 0;

                    switch numel(currSensorAttachmentArray)
                        case 3
                            P = P_data{i}(j, :);
                            P_use = ~isempty(find(P_ind{i} == j, 1));
                            useEntry = [norm(P)] ~= 0;
                            keepEntry = [P_use] ~= 0;
                            if sum(useEntry) == 1 && sum(keepEntry) == 1
                                R_0_imu_curr = eye(3);
                                t_0_imu_curr = mean([P], 1);
                                
                                T = eye(4);
                                T(1:3, 1:3) = R_0_imu_curr;
                                T(1:3, 4) = t_0_imu_curr;
                                
                                T_f_imu_raw{j} = T;
                                T_f_0_avg{j} = T_f_0;
                                R_f_imu{j} = T_f_0*T;
                            else
                                % leave as default []
                            end

                        case {5, 6}
                            P = P_data{i}(j, :);
                            Q = Q_data{i}(j, :);
                            R = R_data{i}(j, :);

                            P_use = ~isempty(find(P_ind{i} == j, 1));
                            Q_use = ~isempty(find(Q_ind{i} == j, 1));
                            R_use = ~isempty(find(R_ind{i} == j, 1));

                            useEntry = [norm(P) norm(Q) norm(R)] ~= 0;
                            keepEntry = [P_use Q_use R_use] ~= 0;

                            if sum(useEntry) >= 3 && sum(keepEntry) >= 3
                                %Rotation
                                [R_0_imu_curr, ~] = points2rot(P,Q,R);
                                
                                %Translation
                                posMarkers = [P; Q; R];
                                t_0_imu_curr = mean(posMarkers, 1);
%                                 [B, TF] = rmoutliers(posMarkers);
                                
                                distPQ = abs(diff(normVector([P; Q])));
                                distPR = abs(diff(normVector([P; R])));
                                distQR = abs(diff(normVector([Q; R])));
                                distAll = [distPQ; distPR; distQR];
%                                  [B, TF] = rmoutliers(distAll);
                                
                                distThres = 0.10; % unlikely that any of the clusters are more than 10 cm away from each other
                                if sum(distAll > distThres)
                                    % leave as default [], and don't use
                                   clearCurrInd = 1;
                                    
                                else
                                    % combine into transformation matrix
                                    T = eye(4);
                                    T(1:3, 1:3) = R_0_imu_curr;
                                    T(1:3, 4) = t_0_imu_curr;
                                    
                                    T_f_imu_raw{j} = T;
                                    T_f_0_avg{j} = T_f_0;
                                    R_f_imu{j} = T_f_0*T;
                                end
                            else
                                % leave as default []
                               clearCurrInd = 1;
                            end
                            
                        otherwise
                            T_f_imu_raw{1} = [];
                            T_f_0_avg{1} = [];
                            R_f_imu{1} = [];
                    end
                    
                    if clearCurrInd
                        for k = 3:length(currSensorAttachmentArray)
%                             indFlatSensor = find(ismember(flatSensorSecondaryLastCell, currSensorAttachmentArray{k}));
                            sensorLengthUseCell{i}(j) = 0;
                        end
                    end
                end
               
                % SVD the rotations together
                T_f_raw_svdAvg = meanSVD_T(T_f_imu_raw);
                T_f_0_svdAvg = meanSVD_T(T_f_0_avg);
                T_f_imu_svdAvg = meanSVD_T(R_f_imu);
                
                if sum(addMarker(i, :)) == expectedMarkers(i) && ~isempty(T_f_raw_svdAvg)
                    transform(i) = imuTrans(currSensorAttachmentArray{2}, frame, ...
                        T_f_raw_svdAvg, T_f_0_svdAvg, T_f_imu_svdAvg);
                else
                    transform(i) = imuTrans(currSensorAttachmentArray{2}, frame);
                end
            end
            
            for i = 1:length(sensorLengthUseCell)
                sensorLengthUseCell{i} = sensorLengthUseCell{i}(sensorLengthUseCell{i} > 0); % remove all the zeros that got left by clearCurrInd
            end
        end
        
        function addJointConstraints(obj, ekf)            
            if isempty(obj.jointRangeOfMotion)
                return
            end
            
            allFrameNames = {obj.model.joints.name};
            for i = 1:length(obj.jointRangeOfMotion)
                indFrame = find(ismember(allFrameNames, obj.jointRangeOfMotion{i}{1}) == 1); 
                
                if ~isempty(indFrame)
                    lb = deg2rad(obj.jointRangeOfMotion{i}{2}(1));
                    ub = deg2rad(obj.jointRangeOfMotion{i}{2}(2));
                    ekf.upper_state_bound(indFrame) = ub + deg2rad(obj.jointRangeOfMotionTolerance); % add margin of error
                    ekf.lower_state_bound(indFrame) = lb - deg2rad(obj.jointRangeOfMotionTolerance);
                else
                    fprintf('rlModelInstance: Joint range available for %s but not found in model\n', obj.jointRangeOfMotion{i}{1});
                end
            end
        end
        
        function markers = collateSensorMarkers(obj) 
            markers = {};
            for i = 1:size(obj.sensorAttachmentArray, 1)
                currSensorAttachmentArray = obj.sensorAttachmentArray{i};
                for j = 3:size(currSensorAttachmentArray, 2)
                    markers = [markers currSensorAttachmentArray(j)];
                end
            end
        end
        
        function loadDemographics(obj)
            for i = 1:length(obj.subjectDemographicArray)
                if obj.subjectDemographicArray{i}{1} == obj.subjectId
                    obj.body_height = obj.subjectDemographicArray{i}{2};
                    obj.body_weight = obj.subjectDemographicArray{i}{3};
                    obj.gender = obj.subjectDemographicArray{i}{4};
                end
            end
            
            if isempty(obj.body_height)
                obj.body_height = 1.70;
                obj.body_weight = 60;
                obj.gender = 'm';
            end
        end
        
        function ind = findMarkerInd(obj, markerNames) 
            % given array of marker names, find the corresponding ind in
            % markerPos
            ind = [];
            for i = 1:length(markerNames)
                newInd = find(strcmpi(obj.measurement_label, markerNames{i}));
                ind = [ind; newInd];
            end
        end
        
         function [kinematicTransform, dynamicTransform] = assembleKinDynTransform(obj, linkFrameStr, sourceFrameStr, t, m, com, I)
            kinematicTransform.frameName = linkFrameStr;
             if ~isempty(t)
                kinematicTransform.t = t;
             else
                kinematicTransform.t = [];
             end
            
            dynamicTransform.frameName = sourceFrameStr;
            if ~isempty(m)
                dynamicTransform.m = m;
                dynamicTransform.com = com;
                dynamicTransform.I   = I;
            else
                dynamicTransform.m = [];
                dynamicTransform.com = [];
                dynamicTransform.I   = [];
            end
         end
        
         function sensorTransform = assembleSenTransform(obj, sensorMarkerStr, attachmentFrameStr, decorators, T)
             sensorTransform.frameName = attachmentFrameStr;
             sensorTransform.sensorName = sensorMarkerStr;
             sensorTransform.decorators = decorators;
             sensorTransform.t = T;
         end
         
         function sensorTransform = assembleSenTransform3dTrc(obj, sensorMarkerStr, attachmentFrameStr, decorators, R, t)
             sensorTransform.frameName = attachmentFrameStr;
             sensorTransform.sensorName = sensorMarkerStr;
             sensorTransform.decorators = decorators;
             
             T = eye(4);
             T(1:3, 1:3) = R;
             T(1:3, 4) = t;
             sensorTransform.t = T;
         end
         
         function sensorTransform = assembleSensorTransformEmpty(obj, sensorMarkerStr, attachmentFrameStr, decorators)
             sensorTransform.frameName = attachmentFrameStr;
             sensorTransform.sensorName = sensorMarkerStr;
             sensorTransform.t = [];
             sensorTransform.decorators = decorators;
         end
        
        function addKinDynTransformToModel(obj, mdl, kinematicTransform, dynamicTransform)
            for i = 1:length(kinematicTransform)
                if ~isempty(kinematicTransform(i).t)
                    obj.setKinematicTransform(mdl, kinematicTransform(i));
                end
            end
            
             for i = 1:length(dynamicTransform)
                if ~isempty(dynamicTransform(i).m)
                    obj.setDynamicTransform(mdl, dynamicTransform(i));
                end
             end
        end
        
        function addSenTransformToModel(obj, mdl, sensorTransform, sensorList)
            indSensors = 0; 
            
            for i = 1:length(sensorTransform)
                if nargin == 4
                    % check to see if the incoming sensor is part of the list
                    % set by the model
                    indFrame = find(ismember(sensorList, sensorTransform(i).sensorName) == 1);
                    
                    if isempty(indFrame)
                        continue;
                    end
                end
                
                indSensors = indSensors + 1;
%                 obj.sensorNameFull{indSensors} = sensorTransform(i).sensorName;
                
                fprintf('Adding sensor %s (%u) to model...\n', sensorTransform(i).sensorName, i);
                
                if ~isempty(sensorTransform(i).t)
                    obj.addSensorTransform(mdl, sensorTransform(i));
%                     obj.sensorAddStatus(indSensors) = 1;
                else
%                     obj.sensorAddStatus(indSensors) = 0;
                end
            end
        end
        
        function keepInds = addSecondaryTransformToModel(obj, mdl, sensorSecondaryTransform, sensorList)
%             lenSen = length(sensorTransform);
            lenSen = length(mdl.sensors);
            indSensors = 0;
            keepInds = [];
            for i = 1:length(sensorSecondaryTransform)
                 if nargin == 4
                    % check to see if the incoming sensor is part of the list
                    % set by the model
                    indFrame = find(ismember(sensorList, sensorSecondaryTransform(i).sensorName) == 1);
                    
                    if isempty(indFrame)
                        continue;
                    end
                 end
                
                 indSensors = indSensors + 1;
                 keepInds = [keepInds i];
%                  obj.sensorNameFull{lenSen + indSensors} = sensorSecondaryTransform(i).sensorName;
                 
                 if ~isempty(sensorSecondaryTransform(i).t)
                     obj.addSensorTransform(mdl, sensorSecondaryTransform(i));
%                      obj.sensorAddStatus(lenSen + indSensors) = 1;
                 else
%                      obj.sensorAddStatus(lenSen + indSensors) = 0;
                 end
            end
        end
        
        function setKinematicTransform(obj, mdl, kinematicTransform)
            allFrameNames = {mdl.transforms.name};
            indFrame = find(ismember(allFrameNames, kinematicTransform.frameName) == 1);
            
            if ~isempty(indFrame)
%                 trans = mdl.transforms(indFrame);
%                 f_in = trans.frame_in;
%                 f_in_rot =  f_in.t(1:3,1:3)';
%                 cartOffset = kinematicTransform.t(1:3, 4);
%                 trans.t(1:3,4) = f_in_rot * cartOffset;
                
                mdl.transforms(indFrame).t = kinematicTransform.t;
            else
                fprintf('WARNING: setKinematicTransform did not find frame: %s\n', kinematicTransform.frameName);
            end
        end
        
        function setDynamicTransform(obj, mdl, dynamicTransform)
            allFrameNames = {mdl.bodies.name};
            indFrame = find(ismember(allFrameNames, dynamicTransform.frameName) == 1);
            
            if ~isempty(indFrame)
                mdl.bodies(indFrame).m   = dynamicTransform.m;
                mdl.bodies(indFrame).com = dynamicTransform.com;
                mdl.bodies(indFrame).I   = dynamicTransform.I;
            else
                fprintf('WARNING: setDynamicTransform did not find frame: %s\n', dynamicTransform.frameName);
            end
        end
        
        function addSensorTransform(obj, mdl, sensorTransform)
            % check to make sure the frame exist
            f = mdl.getFrameByName(sensorTransform.frameName);
            
            if isempty(f)
                fprintf('WARNING: addSensorTransform did not find frame: %s\n', sensorTransform.frameName);
            end
            
            M1 = SensorCore(sensorTransform.sensorName);
            for i = 1:length(sensorTransform.decorators)
                M1.addDecorator(sensorTransform.decorators{i});
            end
            mdl.addSensor(M1, sensorTransform.frameName, sensorTransform.t);
            
            clear M1;
        end
        
%         function addSensorTransformEkfModel(obj, sensorTransform)
%             M1 = SensorCore(sensorTransform.sensorName);
%             for i = 1:length(sensorTransform.decorators)
%                 M1.addDecorator(sensorTransform.decorators{i});
%             end
%             obj.model_ekf.addSensor(M1, sensorTransform.frameName, sensorTransform.t);
%             
%             clear M1;
%         end
        
        function checkSymmetry(obj)
            % if the mode is X00, then check these sides for symmetry
            % use the actual function to get the signage
        end
        
        function sensorNameFull = getSensorNameFull(obj)
            if ~isempty(obj.model_ekf) && length(obj.model_ekf.sensors) > 0
                sensorNameFull = {obj.model_ekf.sensors.name};
            else
                sensorNameFull = [];
            end
        end
    end
    
    methods(Static = true)
        function [sourceFrameStr, targetFrameStr, sourceMarkerStr, targetMarkerStr] = linkLengthArrayString(currLinkLengthArrangement)
            targetFrameStr = currLinkLengthArrangement{1};
            sourceFrameStr = currLinkLengthArrangement{2};
            sourceMarkerStr = currLinkLengthArrangement{3};
            targetMarkerStr = currLinkLengthArrangement{4};
        end
        
        function [sourceFrameStr, sourceMarkerStr, targetMarkerStr] = sensorLengthArrayString(currLinkLengthArrangement)
            sourceFrameStr = currLinkLengthArrangement{1};
            sourceMarkerStr = currLinkLengthArrangement{2};
            targetMarkerStr = currLinkLengthArrangement{3};
        end
        
        function [sourceMarkerData, targetMarkerData] = markerLengthData(sourceMarkerStr, targetMarkerStr, dataInstance)
            % calculate link offsets
            sourceMarkerData = dataInstance.data.(sourceMarkerStr);
            targetMarkerData = dataInstance.data.(targetMarkerStr);
        end
        
        function indsMarkersToKeep = checkMarkerNorm(markers)
            normMarkers = normVector(markers);
            indsMarkersToKeep = find(normMarkers > 0);
        end
        
        function [sourceFrameStr, targetFrameStr, linkLength, linkVectorMean, linkVectorStd] = linkLengthArrayData(currLinkLengthArrangement, linkLengthInd, dataInstance)
            [sourceFrameStr, targetFrameStr, sourceMarkerStr, targetMarkerStr] = rlModelInstance.linkLengthArrayString(currLinkLengthArrangement);
            [linkLength, linkVectorMean, linkVectorStd] = rlModelInstance.combineMarkers(sourceMarkerStr, targetMarkerStr, linkLengthInd, dataInstance);
        end
        
        function [sourceFrameStr, targetMarkerStr, linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData, useIntersect] = sensorLengthArrayData(currLinkLengthArrangement, linkLengthInd, dataInstance)
            [sourceFrameStr, sourceMarkerStr, targetMarkerStr] = rlModelInstance.sensorLengthArrayString(currLinkLengthArrangement);
            [linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData, useIntersect] = rlModelInstance.combineMarkers(sourceMarkerStr, targetMarkerStr, linkLengthInd, dataInstance);
        end
        
        function [linkLength, linkVectorMean, linkVectorStd, sourceMarkerData, targetMarkerData, useIntersect2] = combineMarkers(sourceMarkerStr, targetMarkerStr, linkLengthInd, dataInstance)
            [sourceMarkerData, targetMarkerData] = rlModelInstance.markerLengthData(sourceMarkerStr, targetMarkerStr, dataInstance);
            
            % check the markers outlined in linkLengthInd for empty data
            sourceIndsNormKeep = rlModelInstance.checkMarkerNorm(sourceMarkerData);
            targetIndsNormKeep = rlModelInstance.checkMarkerNorm(targetMarkerData);
            
            % figure out the intersection between markers that are not at zero
            % and the indUse array being passed in
            useIntersect1 = intersect(sourceIndsNormKeep, targetIndsNormKeep);
            useIntersect2 = intersect(useIntersect1, linkLengthInd);
            
            % average together all the marker data
            linkLengthArray = targetMarkerData(useIntersect2, :) - sourceMarkerData(useIntersect2, :);
%             linkVectorMean = mean(linkLengthArray, 1);
            
            normLen = normVector(linkLengthArray);
            linkLength = mean(normLen, 1);
            
            if isempty(normLen)
                linkVectorMean = NaN*ones(3, 1);
            elseif normLen(1) > 0
                scaleNorm = normLen(1)/mean(normLen);
                linkVectorMean = linkLengthArray(1, :)/scaleNorm;
            else
                linkVectorMean = mean(linkLengthArray(:, :), 1);
            end
            
            linkVectorStd = std(linkLengthArray, 1);
        end
              
        function sensorAttachment = flattenSensorAttachmentArray(attachmentArray)
            indSensors = 0;
            sensorAttachment = {};
            for i = 1:size(attachmentArray, 1)
                currSensorAttachmentArray = attachmentArray{i};
                
                for j = 3:size(currSensorAttachmentArray, 2)
                    indSensors = indSensors + 1;
                    sensorAttachment{indSensors, 1} = [currSensorAttachmentArray(1) currSensorAttachmentArray(2) currSensorAttachmentArray(j)];
                end
            end
        end
    end
end