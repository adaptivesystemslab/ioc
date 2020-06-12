classdef rlFeatureSet < handle % < rlModel
    properties
        % subject specific data
        subjectParam = [];
        algorithmParam = [];
        
        model;
        dt;
        
        % performing inverse kinematics
        measurement_labels = []; % this may be alphabetical, or in sensor order
        measurement_input = [];
        measurement_output = [];
        measurement_mask = [];
        measurement_ekf_match = [];
        measurement_input_match = [];
        measurement_output_match = [];

        % estimates
        ekf = [];
        ekf_states = [];
        ekf_eventhandler = [];
        ekf_match = [];
        ekf_matchLabels = []; % expect the input to be alphabetical
        segments = [];
        
        % other framing information
        baseFrame;
        frameData;
        
        visSkipRate = 1; 
    end
    
    methods
        function obj = rlFeatureSet()
            
        end
        
        function setModelParam(obj, mdl, dt, algorithmParam)
            obj.model = mdl;
            obj.dt = dt;
            
            obj.algorithmParam = algorithmParam;  
%             if obj.algorithmParam.ekfRun(end) > numel(pModelInstance.time)
%                 obj.algorithmParam.ekfRun = obj.algorithmParam.ekfRun(1):numel(pModelInstance.time);
%             end       
        end
        
        function setEKFParams(obj, rlModelInstance, ekfType, tuningParams)
            % EKF model            
            obj.ekf = setupEKF(obj.model, obj.dt, ekfType, tuningParams);
            obj.joint_labels = {obj.model.joints.name};
            
            ekf_eventhandler = EventHandler();
            obj.ekf.addlistener('matched_callback',@ekf_eventhandler.matched_marker);
            obj.ekf.addlistener('missing_callback',@ekf_eventhandler.missing_marker);
            obj.ekf.addlistener('swap_callback',@ekf_eventhandler.swap_marker);
            obj.ekf_eventhandler = ekf_eventhandler;
            
            % set up constraints
            if ~isempty(rlModelInstance)
                rlModelInstance.addJointConstraints(obj.ekf);
            end
            
            % update marker swapping parameter
            obj.ekf.update_smt_matrix(1, tuningParams.markerSwapTolerance); % marker position
            obj.ekf.update_smt_matrix(17, tuningParams.markerSwapTolerance); % marker position/velocity
        end
        
        function measurements = calcMeasurementFromState(obj, featureSet, dataInstance)
            mdl = obj.model;
            
            numTimeLength = size(featureSet.q, 1);
            obj.measurement_output = SensorMeasurement(numTimeLength, length(mdl.sensors));
            obj.measurement_labels = {mdl.sensors.name};
            
            obj.q = featureSet.q;
            obj.dq = featureSet.dq;
            obj.ddq = featureSet.ddq;
            obj.joint_labels = featureSet.joint_labels;
            
            indArray = 1:size(obj.q, 1);
            
            for j = 1:length(mdl.sensors)
                sensorArray = SensorMeasurement(numTimeLength, 1);
                [sensorArray.type] = deal(mdl.sensors(1).binary_type);
                [sensorArray.size] = deal(numel(mdl.sensors(1).measurement));
                obj.measurement_output(:, j) = sensorArray;
            end
            
            % visualize the reconstruction
            if obj.algorithmParam.visualize
                vis = rlVisualizer('vis',640,480);
                mdl.forwardPosition();
                vis.addModel(mdl);
                vis.update();
            end
            
            f = rlFrame.empty();
            for j = 1:length(obj.frameData)
                f(end+1) = mdl.getFrameByName(obj.frameData(j).name);
            end
            
            % run EKF
            for i = 1:numTimeLength
                mdl.position = obj.q(i, :);
                mdl.velocity = obj.dq(i, :);
                mdl.acceleration = obj.ddq(i, :);

                mdl.forwardPosition();
                mdl.forwardVelocity();
                mdl.forwardAcceleration();
    
                if (i == 1 || mod(i, obj.visSkipRate) == 0)
                    fprintf('Frame %u of %u\n', i, numTimeLength);
                    
                    if obj.algorithmParam.visualize 
                        if ~isempty(dataInstance)
                            % visualize markers
                            applyMarkersToVisualization(vis, mdl, dataInstance, indArray(i));
                        end
                        
                        vis.update();
                        pause(0.01);
                    end
                end
                
                 if ~isempty(obj.frameData)
                    mdl.forwardPosition();
                    mdl.forwardVelocity();
                    mdl.forwardAcceleration();
                    
                    for j = 1:length(obj.frameData)
                        obj.frameData(j).position(i, :) = reshape(f(j).t, 16, 1);
                        obj.frameData(j).velocity(i, :) = f(j).v;
                        obj.frameData(j).acceleration(i, :) = f(j).a;
                    end
                end
                
                obj.measurement_output(i, :) = SensorMeasurement(mdl.sensors);
            end
            
            measurements = obj.measurement_output;
        end
        
        function calcFkFromState(obj, featureSet, dataInstance)
            mdl = obj.model;
            
            numTimeLength = size(featureSet.q, 1);
            obj.measurement_output = SensorMeasurement(numTimeLength, length(mdl.sensors));
            obj.measurement_labels = {mdl.sensors.name};
            
            obj.q = featureSet.q;
            obj.dq = featureSet.dq;
            obj.ddq = featureSet.ddq;
            obj.joint_labels = featureSet.joint_labels;
            
            indArray = 1:size(obj.q, 1);
            
            if isempty(obj.algorithmParam) || isempty(obj.algorithmParam.externalBaseFrameData)
                externalBaseFrameData = [];
            else
                sourceFrameInjectionStr = 'body_base';
                allNamesExternalData = {obj.algorithmParam.externalBaseFrameData.name};
                indFrameExternalData = find(ismember(allNamesExternalData, sourceFrameInjectionStr) == 1);
                
                targetFrameInjectionStr = 'frame_6dof_revx2_to_body_base';
                allNamesTransforms = {mdl.transforms.name};
                indFrameTransforms = find(ismember(allNamesTransforms, targetFrameInjectionStr) == 1);
                
                if ~isempty(indFrameExternalData)
                    externalBaseFrameData.position = obj.algorithmParam.externalBaseFrameData(indFrameExternalData).position;
                    externalBaseFrameData.velocity = obj.algorithmParam.externalBaseFrameData(indFrameExternalData).velocity;
                    externalBaseFrameData.acceleration = obj.algorithmParam.externalBaseFrameData(indFrameExternalData).acceleration;
                else
                    externalBaseFrameData = [];
                end
            end
			
            % visualize the reconstruction
            if obj.algorithmParam.visualize
                vis = rlVisualizer('vis',640,480);
                mdl.forwardPosition();
                vis.addModel(mdl);
                vis.update();
            end
            
            % run EKF
            for i = 1:numTimeLength
                mdl.position = obj.q(i, :);
                mdl.velocity = obj.dq(i, :);
                mdl.acceleration = obj.ddq(i, :);
                
                if ~isempty(externalBaseFrameData)
                    mdl.transforms(indFrameTransforms).t = reshape(externalBaseFrameData.position(i, :), 4, 4);
                    mdl.forwardPosition();
                end
                
                mdl.forwardPosition();
                mdl.forwardVelocity();
                mdl.forwardAcceleration();
                
                if (i == 1 || mod(i, obj.visSkipRate) == 0)
                    fprintf('Frame %u of %u\n', i, numTimeLength);
                    
                    if obj.algorithmParam.visualize
                        if ~isempty(dataInstance)
                            % visualize markers
                            applyMarkersToVisualization(vis, mdl, dataInstance, indArray(i));
                        end
                        
                        vis.update();
                        pause(0.01);
                    end
                end
            end
        end
        
%         % calculate the EKF states from EKF 
%         function ekf_states = calcStatesFromTrc(obj, measurement_input, initPos, ...
%                 dataInstance, modelInstance, algorithmParam, indArray, labelArray, matchMatrix)
%             mdl = obj.model;
%             obj.joint_labels = mdl.joints.names;
%             
%             if isempty(initPos)
%                 initPos = zeros(length(mdl.position), 1);
%             end
%             
%             if isempty(indArray)
%                 indArray = 1:size(measurement_input, 1);
%             end
%             
%              if exist('matchMatrix', 'var') && ~isempty(matchMatrix)
%                 matches = matchMatrix(1:length(mdl.sensors), :);
%             else
%                 matches = [1:length(mdl.sensors); 1:length(mdl.sensors)]';
%             end
%             
%             obj.measurement_input = measurement_input;
%             obj.measurement_output = measurement_input(:, 1:length(mdl.sensors)); % will be replaced later in the fct
%             obj.measurement_labels = {mdl.sensors.name};
%             obj.ekf_matchLabels = labelArray;
%             
%             % initialize the array
%             numTimeLength = size(measurement_input, 1);
%             obj.ekf_states = zeros(numTimeLength, numel(obj.ekf.state));
%             obj.ekf.state(1:length(mdl.position)) = initPos(1:length(mdl.position)); % for the first timestep, set it to the current model value
%             obj.ekf.state(length(mdl.position)+1:end) = 0;
%             
%             mdl.position = initPos;
%             mdl.forwardPosition();
%             
%             % visualize the reconstruction
%             if obj.algorithmParam.visualize
%                 vis = rlVisualizer('vis',640,480);
%                 vis.addModel(mdl);
%                 vis.update();
%             end
% 
%             % run EKF
%             for i = 1:numTimeLength
%                 mes_obj = measurement_input(i, :);
%                 u = obj.dt;
%    
%                 if obj.algorithmParam.ekfForceMatch == 1
%                     cur_ass = obj.ekf.run_iteration(u, mes_obj, matches);   % calculate the IK
%                 elseif obj.algorithmParam.ekfForceMatch == -1 % match only the first frame
%                     cur_ass = obj.ekf.run_iteration(u, mes_obj, matches);   % calculate the IK
%                     obj.algorithmParam.ekfForceMatch = 0;
%                 else
%                     cur_ass = obj.ekf.run_iteration(u, mes_obj);
%                 end
%                 
%                 if (i == 1 || mod(i, obj.visSkipRate) == 0)
%                     fprintf('Frame %u of %u\n', i, numTimeLength);
%                     
%                     if obj.algorithmParam.visualize 
%                         if ~isempty(dataInstance)
%                             % visualize markers
%                             applyMarkersToVisualization(vis, mdl, dataInstance, indArray(i), cur_ass.mes, matchMatrix);
%                         end
%                         
%                         vis.update();
%                         pause(0.01);
%                     end
%                 end
%                 
%                 if i == 1
%                     if 0
%                         rlFeatureSet.listMarkerMatches(matches(:, 2)');
%                         rlFeatureSet.listMarkerMatches(cur_ass.getMesArray);
%                     end
%                     
% %                     fprintf('Match before iteration: \n');
% %                     rlFeatureSet.listMarkerMatches(mdl, labelArray, matches(:, 2)');
%                     fprintf('Match after iteration: \n');
%                     rlFeatureSet.listMarkerMatches(mdl, labelArray, cur_ass.getMesArray);
%                 end
% 
%                 obj.ekf_match(i, :) = cur_ass.mes';
%                 obj.ekf_states(i, :) = obj.ekf.state; % pull out the IK
%                 obj.measurement_output(i, :) = obj.ekf.z_predict;     % pull out the FK
%             end
%             
%             ekf_states = obj.ekf_states;
%             
%             if obj.algorithmParam.visualize
%                 clear vis;
%             end
%         end
        
          % calculate the EKF states from EKF 
        function ekf_states = calcStatesFromTrc_BaseFrame(obj, measurement_input, qInit, dqInit, ddqInit, ...
                dataInstance, modelInstance, indArray, labelArray, matchMatrix)
            mdl = obj.model;
            obj.joint_labels = {mdl.joints.name};
            
            if isempty(qInit)
                qInit = zeros(length(mdl.position), 1);
                dqInit = zeros(length(mdl.velocity), 1);
                ddqInit = zeros(length(mdl.acceleration), 1);
            end
            
            if isempty(dqInit)
                dqInit = zeros(length(mdl.velocity), 1);
                ddqInit = zeros(length(mdl.acceleration), 1);
            end
            
            if isempty(indArray)
                indArray = 1:size(measurement_input, 1);
            end
            
            if exist('matchMatrix', 'var') && ~isempty(matchMatrix)
                inds = find(matchMatrix(:, 1) > 0);
                matches = matchMatrix(inds, :);
            else
                matches = [1:length(mdl.sensors); 1:length(mdl.sensors)]';
            end
            
            if isempty(obj.algorithmParam) || isempty(obj.algorithmParam.externalBaseFrameData)
                externalBaseFrameData = [];
            else
                sourceFrameInjectionStr = 'body_base';
                allNamesExternalData = {obj.algorithmParam.externalBaseFrameData.name};
                indFrameExternalData = find(ismember(allNamesExternalData, sourceFrameInjectionStr) == 1);
                
                targetFrameInjectionStr = 'frame_6dof_revx2_to_body_base';
                allNamesTransforms = {mdl.transforms.name};
                indFrameTransforms = find(ismember(allNamesTransforms, targetFrameInjectionStr) == 1);
                
                if ~isempty(indFrameExternalData)
                    externalBaseFrameData.position = obj.algorithmParam.externalBaseFrameData(indFrameExternalData).position;
                    externalBaseFrameData.velocity = obj.algorithmParam.externalBaseFrameData(indFrameExternalData).velocity;
                    externalBaseFrameData.acceleration = obj.algorithmParam.externalBaseFrameData(indFrameExternalData).acceleration;
                else
                    externalBaseFrameData = [];
                end
            end
                 
            % setup preallocation for the output. make sure the output has the
            % same size as the rmmodelinstance specifications, not with the
            % mdl.sensor length, just in case some sensors were not added to the
            % model since the initialization data is poor
            obj.measurement_input = measurement_input;
            obj.measurement_output = measurement_input(:, 1:length(mdl.sensors)); % will be replaced later in the fct
            obj.measurement_labels = {mdl.sensors.name};
            obj.ekf_matchLabels = labelArray;
            obj.ekf_match = [];
            
            % initialize ekf and its storage array
            numTimeLength = size(measurement_input, 1);
            statePosInd = [1:length(mdl.position)] + length(mdl.position)*0;
            stateVelInd = [1:length(mdl.position)] + length(mdl.position)*1;
            stateAccInd = [1:length(mdl.position)] + length(mdl.position)*2;
            obj.ekf_states = zeros(numTimeLength, numel(obj.ekf.state));
            obj.ekf.state(statePosInd) = qInit; % for the first timestep, set it to the current model value
            if length(obj.ekf.state) >= max(stateVelInd)
                obj.ekf.state(stateVelInd) = dqInit;
            end
            if length(obj.ekf.state) >= max(stateAccInd)
                obj.ekf.state(stateAccInd) = ddqInit;
            end
            
            % initialize the model in a similar fashion
            mdl.position = qInit;
            mdl.velocity = dqInit;
            mdl.acceleration = ddqInit;
            mdl.forwardPosition();
            mdl.forwardVelocity();
            mdl.forwardAcceleration();
            
            % visualize the reconstruction
            if obj.algorithmParam.visualize
                vis = rlVisualizer('vis',640,480);
                vis.addModel(mdl);
                applyMarkersToVisualization(vis, mdl, dataInstance, 1, [], matchMatrix);
                vis.update();
            end
            
            f = rlFrame.empty();
            for j = 1:length(obj.frameData)
                fprintf('Locating frame %s\n', obj.frameData(j).name);
                f(j) = mdl.getFrameByName(obj.frameData(j).name);
                f(j); % check to see if the connection has been made
            end
             
            % run EKF
            for i = 1:numTimeLength
                mes_obj = measurement_input(i, :);
                u = obj.dt;
                
                if ~isempty(externalBaseFrameData)
                    mdl.transforms(indFrameTransforms).t = reshape(externalBaseFrameData.position(i, :), 4, 4);
                    mdl.forwardPosition();
                end
                
                switch obj.algorithmParam.ekfForceMatch
                    case 1 
                        % always force a match
                        [cur_ass, cur_ass_out] = obj.ekf.run_iteration(u, mes_obj, matches);   % calculate the IK

                    case -1 
                        % force match only in first frame
                        [cur_ass, cur_ass_out]  = obj.ekf.run_iteration(u, mes_obj, matches);   % calculate the IK
                        obj.algorithmParam.ekfForceMatch = 0;

                    case -2
                        % force a match only if the marker is in a reasonable
                        % location, ie not at 10000
                        matchesMod = matches;
                        for j = 1:length(mes_obj)
                            currMesCheck = mes_obj(j).getMesArray;
                            if abs(currMesCheck(1)) > 1e4
                                findMatchInd = find(matchesMod(:, 2) == j);
                                if ~isempty(findMatchInd)
                                    matchesMod = [matchesMod(1:(findMatchInd-1), :); matchesMod(findMatchInd+1:end, :)];
                                end
                            end
                        end
                        [cur_ass, cur_ass_out]  = obj.ekf.run_iteration(u, mes_obj, matchesMod);   % calculate the IK

                    otherwise
                        % let ekf determine the match
                        [cur_ass, cur_ass_out]  = obj.ekf.run_iteration(u, mes_obj);
                end
                
                % flatten cur_ass if need to. just keep the first one for
                % now
                if length(cur_ass_out) > 1
                    cur_ass = cur_ass_out(1);
                end
                
                if (i == 1 || mod(i, obj.visSkipRate) == 0)
%                     fprintf('Frame %u of %u\n', i, numTimeLength);
                    
                    if obj.algorithmParam.visualize 
                        if ~isempty(dataInstance)
                            % visualize markers
                            applyMarkersToVisualization(vis, mdl, dataInstance, indArray(i), cur_ass.mes, matchMatrix);
                        end
                        
                        vis.update();
                        pause(0.01); % break
                    end
                end
                
                if i == 1
                    if 0
                        rlFeatureSet.listMarkerMatches(mdl, labelArray, matches(:, 2)');
                        rlFeatureSet.listMarkerMatches(mdl, labelArray, cur_ass.getMesArray);
                    end
                    
%                     fprintf('Match before iteration: \n');
%                     rlFeatureSet.listMarkerMatches(mdl, labelArray, matches(:, 2)');
                    fprintf('Match after first iteration: \n');
                    rlFeatureSet.listMarkerMatches(mdl, labelArray, cur_ass.getMesArray);
                end
                
                if ~isempty(obj.frameData)
                    mdl.forwardPosition();
                    mdl.forwardVelocity();
                    mdl.forwardAcceleration();
                    
                    for j = 1:length(obj.frameData)
                        obj.frameData(j).position(i, :) = reshape(f(j).t, 16, 1);
                        obj.frameData(j).velocity(i, :) = f(j).v;
                        obj.frameData(j).acceleration(i, :) = f(j).a;
                    end
                end
                
                obj.ekf_match(i, :) = cur_ass.mes';
                obj.ekf_states(i, :) = obj.ekf.state; % pull out the IK
                obj.measurement_output(i, :) = obj.ekf.z_predict;     % pull out the FK
            end
            
            ekf_states = obj.ekf_states;
            
            if 0
                vis = rlVisualizer('vis',640,480);
                vis.addModel(mdl);
                applyMarkersToVisualization(vis, mdl, dataInstance, indArray(i), cur_ass.mes, matchMatrix);
                vis.update();
            end
            
            if exist('vis', 'var')
                clear vis;
            end
        end
  
         function ekf_states = calcStatesFromTrc_BaseFrame3(obj, measurement_input, initPose_R, initPose_L, ...
                dataInstance, modelInstance, algorithmParam, indArray, labelArray, matches, mdl_R, mdl_L, algorithmParam_R, algorithmParam_L)
%             mdl = obj.model;
            
            modelPosition = [mdl_R.position; mdl_L.position];
            modelJoints = {};
            for i = 1:length(mdl_R.joints)
                modelJoints{i} = [mdl_R.joints(i).name '_mdl_R'];
                modelJoints{i+length(mdl_R.joints)} = [mdl_L.joints(i).name '_mdl_L'];
            end
%             modelJoints = [{mdl_R.joints.name} {mdl_L.joints.name}];
            modelSensors = [mdl_R.sensors([mdl_R.sensors.binary_type] == 6) mdl_L.sensors([mdl_L.sensors.binary_type] == 6)];
            
            % setup preallocation for the output. make sure the output has the
            % same size as the rmmodelinstance specifications, not with the
            % mdl.sensor length, just in case some sensors were not added to the
            % model since the initialization data is poor
            obj.measurement_input = measurement_input;
            obj.measurement_output = measurement_input(:, 1:length(modelSensors)); % will be replaced later in the fct
            obj.measurement_labels = {modelSensors.name};
            obj.ekf_matchLabels = labelArray;
            obj.ekf_match = [];
            
            obj.joint_labels = modelJoints;
            
            initPos = [initPose_R initPose_L];
            
            % initialize the array
            numTimeLength = size(measurement_input, 1);
            obj.ekf_states = zeros(numTimeLength, numel(obj.ekf.state));
            obj.ekf.state(1:length(modelPosition)) = initPos(1:length(modelPosition)); % for the first timestep, set it to the current model value
            obj.ekf.state(length(modelPosition)+1:end) = 0;
            
            mdl_R.position = initPose_R;
            mdl_L.position = initPose_L;
            mdl_R.forwardPosition();
            mdl_L.forwardPosition();
            mdl_R.base = algorithmParam_R.baseFrame;
            mdl_L.base = algorithmParam_L.baseFrame;
            
            % visualize the reconstruction
            if obj.algorithmParam.visualize
                vis = rlVisualizer('vis',640,480);
                vis.addModel(mdl_R);
                vis.addModel(mdl_L);
                vis.update();
            end

%     frameNames = {'body_base', ...
%         'frame_rhip_0', 'frame_rhip_6', ...
%         'frame_rknee_0', 'frame_rknee_4', ...
%         'frame_rankle_0', 'frame_rankle_6', ...
%         'frame_lhip_0', 'frame_lhip_6', ...
%         'frame_lknee_0', 'frame_lknee_4', ...
%         'frame_lankle_0', 'frame_lankle_6', ...
%         'frame_lballfoot_end', 'frame_rballfoot_end'};
            f(1) = mdl_R.getFrameByName('body_base');
            f(2) = mdl_R.getFrameByName('frame_rhip_0');
            f(3) = mdl_R.getFrameByName('frame_rhip_6');
            f(4) = mdl_R.getFrameByName('frame_rknee_0');
            f(5) = mdl_R.getFrameByName('frame_rknee_4');
            f(6) = mdl_R.getFrameByName('frame_rankle_0');
            f(7) = mdl_R.getFrameByName('frame_rankle_6');
            f(8) = mdl_L.getFrameByName('frame_lhip_0');
            f(9) = mdl_L.getFrameByName('frame_lhip_6');
            f(10) = mdl_L.getFrameByName('frame_lknee_0');
            f(11) = mdl_L.getFrameByName('frame_lknee_4');
            f(12) = mdl_L.getFrameByName('frame_lankle_0');
            f(13) = mdl_L.getFrameByName('frame_lankle_6');
            f(14) = mdl_L.getFrameByName('frame_lballfoot_end');
            f(15) = mdl_R.getFrameByName('frame_rballfoot_end');
            
            % run EKF
            for i = 1:numTimeLength
                mes_obj = measurement_input(i, :).getMesArray();
                u = obj.dt;
  
                obj.ekf.run_iteration(u, mes_obj');   % calculate the IK
                %         cur_ass = obj.ekf.run_iteration(u, mes_obj);
                
                if (i == 1 || mod(i, obj.visSkipRate) == 0)
                    fprintf('Frame %u of %u\n', i, numTimeLength);
                    
                    if obj.algorithmParam.visualize 
                        if ~isempty(dataInstance)
                            % visualize markers
%                             applyMarkersToVisualization(vis, mdl, dataInstance, indArray(i), cur_ass.mes, matchMatrix);
                            for j = 1:length(mdl_R.sensors)
                                name = ['mdl_R_' mdl_R.sensors(j).name];
                                pos = mdl_R.sensors(j).transform(1:3, 4);
                                vis.addMarker(name, pos, [1 0 0 0]);
                            end
                            for j = 1:length(mdl_L.sensors)
                                name = ['mdl_L_' mdl_L.sensors(j).name];
                                pos = mdl_L.sensors(j).transform(1:3, 4);
                                vis.addMarker(name, pos, [0 1 0 0]);
                            end
                        end
                        
                        vis.update();
                        pause(0.01);
                    end
                end
                
                if ~isempty(obj.frameData)
                    mdl_R.forwardPosition();
                    mdl_R.forwardVelocity();
                    mdl_R.forwardAcceleration();
                    mdl_L.forwardPosition();
                    mdl_L.forwardVelocity();
                    mdl_L.forwardAcceleration();
                    
                    for j = 1:length(obj.frameData)
                        obj.frameData(j).position(i, :) = reshape(f(j).t, 16, 1);
                        obj.frameData(j).velocity(i, :) = f(j).v;
                        obj.frameData(j).acceleration(i, :) = f(j).a;
                    end
                end
                
%                 obj.ekf_match(i, :) = cur_ass.mes';
                obj.ekf_states(i, :) = obj.ekf.state; % pull out the IK
%                 obj.measurement_output(i, :) = obj.ekf.z_predict;     % pull out the FK
            end
            
            ekf_states = obj.ekf_states;
            
            if obj.algorithmParam.visualize
                clear vis;
            end
        end
%           function ekf_states = calcStatesFromTrc_BaseFrame2(obj, measurement_input, initPos, dataInstance, algorithmParam, dataInstance2)
%             mdl = obj.model;
%             
%             if isempty(initPos)
%                 initPos = zeros(length(mdl.position), 1);
%             end
%             
%             if isempty(algorithmParam) || isempty(algorithmParam.externalBaseFrameData)
%                 externalBaseFrameData = [];
%             else
%                 sourceFrameInjectionStr = 'body_world';
%                 targetFrameInjectionStr = 'world';
%                 allNames = {algorithmParam.externalBaseFrameData.name};
%                 
%                 indFrame = find(ismember(allNames, sourceFrameInjectionStr) == 1);
%                 
%                 if ~isempty(indFrame)
%                     externalBaseFrameData.position = algorithmParam.externalBaseFrameData(indFrame).position;
%                     externalBaseFrameData.velocity = algorithmParam.externalBaseFrameData(indFrame).velocity;
%                     externalBaseFrameData.acceleration = algorithmParam.externalBaseFrameData(indFrame).acceleration;
%                 else
%                     externalBaseFrameData = [];
%                 end
%             end
%             
%             obj.measurement_input = measurement_input;
%             obj.measurement_output = measurement_input; % will be replaced later in the fct
%             obj.measurement_labels = {mdl.sensors.name};
%             
%             % initialize the array
%             numTimeLength = size(measurement_input, 1);
%             obj.ekf_states = zeros(numTimeLength, numel(obj.ekf.state));
%             obj.ekf.state(1:length(mdl.position)) = initPos(1:length(mdl.position)); % for the first timestep, set it to the current model value
%             obj.ekf.state(length(mdl.position)+1:end) = 0;
%             
%             % visualize the reconstruction
%             if obj.algorithmParam.visualize
%                 vis = rlVisualizer('vis',640,480);
%                 mdl.forwardPosition();
%                 vis.addModel(mdl);
%                 vis.update();
%             end
%             
%             matches = [1:numel(mdl.sensors)]';
%             matches = [matches matches];
% 
%             f = rlFrame.empty();
%             for j = 1:length(obj.frameData)
%                 f(end+1) = mdl.getFrameByName(obj.frameData(j).name);
%             end
%             
%             % run EKF
%             for i = 1:numTimeLength
%                 mes_obj = measurement_input(i, :);
%                 u = obj.dt;
%                 
% %                 if ~isempty(externalBaseFrameData)
% %                     f = mdl.getFrameByName(targetFrameInjectionStr);
% %                     f.t = reshape(externalBaseFrameData.position(i, :), 4, 4);
% %                     f.v = externalBaseFrameData.velocity(i, :);
% %                     f.a = externalBaseFrameData.acceleration(i, :);
% %                     
% %                     mdl.forwardPosition();
% %                     mdl.forwardVelocity();
% %                     mdl.forwardAcceleration();
% %                 end
%                 
%                 if obj.algorithmParam.ekfForceMatch == 1
%                     cur_ass = obj.ekf.run_iteration(u, mes_obj, matches);   % calculate the IK
%                 elseif obj.algorithmParam.ekfForceMatch == -1 % match only the first frame
%                     cur_ass = obj.ekf.run_iteration(u, mes_obj, matches);   % calculate the IK
%                     obj.algorithmParam.ekfForceMatch = 0;
%                 else
%                     cur_ass = obj.ekf.run_iteration(u, mes_obj);
%                 end
%                 
%                 if (i == 1 || mod(i, obj.visSkipRate) == 0)
%                     fprintf('Frame %u of %u\n', i, numTimeLength);
%                     
%                     if obj.algorithmParam.visualize 
%                         if ~isempty(dataInstance)
%                             % visualize markers
%                             applyMarkersToVisualization(vis, mdl, dataInstance2, i);
%                         end
%                         
%                         % visualize model markers
%                         for j = 1:length(mdl.sensors)
%                             pos = mdl.sensors(j).transform(1:3, 4);
%                             vis.addMarker(['MODEL_' mdl.sensors(j).name], pos, [0.8 0 0 1]);
%                         end
%                         
%                         vis.update();
%                         pause(0.01);
%                     end
%                 end
%                 
%                 if ~isempty(obj.frameData)
%                     mdl.forwardPosition();
%                     mdl.forwardVelocity();
%                     mdl.forwardAcceleration();
%                     
%                     for j = 1:length(obj.frameData)
%                         obj.frameData(j).position(i, :) = reshape(f(j).t, 16, 1);
%                         obj.frameData(j).velocity(i, :) = f(j).v;
%                         obj.frameData(j).acceleration(i, :) = f(j).a;
%                     end
%                    
%                 end
%                 
%                 obj.ekf_match(i, :) = cur_ass.mes';
%                 obj.ekf_states(i, :) = obj.ekf.state; % pull out the IK
%                 obj.measurement_output(i, :) = obj.ekf.z_predict;     % pull out the FK
%             end
%             
%             if obj.algorithmParam.visualize
%                 clear vis;
%             end
%             
%             ekf_states = obj.ekf_states;
%         end
        
        function setFrameSaving(obj, frameNames, numTimeLength)
%             numTimeLength = size(measurement_input, 1);
            
            obj.frameData =  struct('name', frameNames{1}, ...
                'position', zeros(numTimeLength, 16), ...
                'velocity', zeros(numTimeLength, 6), ...
                'acceleration', zeros(numTimeLength, 6));
            for i = 2:length(frameNames)
                obj.frameData(i) = struct('name', frameNames{i}, ...
                'position', zeros(numTimeLength, 16), ...
                'velocity', zeros(numTimeLength, 6), ...
                'acceleration', zeros(numTimeLength, 6));
            end
        end
        
        function saveData(obj, filepath)
            time = obj.time;
            dt = obj.dt;
            ekf_state = obj.ekf_states;
            
            save(filepath, 'time', 'dt', 'ekf_state');
        end
        
        function loadData(obj, filepath)
            loadData = load(filepath);
            
            obj.time = loadData.time;
            obj.dt = loadData.dt;
            obj.ekf_states = loadData.ekf_states;
        end
    end
    
    methods(Static)
        function listMarkerMatches(mdl, labelArray, markerAssignment)
            for j = 1:size(markerAssignment, 2)
                currAss = markerAssignment(j);
                if currAss > 0
                    mdlSenStr = mdl.sensors(j).name;
                    mdlSenNum = j;
                    mesSenStr = labelArray{currAss};
                    mesSenNum = currAss;
                    strMatch = strcmpi(mdlSenStr, mesSenStr);
                else
                    mdlSenStr = mdl.sensors(j).name;
                    mdlSenNum = j;
                    mesSenStr = 'NO MATCH';
                    mesSenNum = 0;
                    strMatch = 0;
                end
                
                fprintf('EKF:MES match (%u) - %s (%u) : %s (%u) \n', ...
                    strMatch, mdlSenStr, mdlSenNum, mesSenStr, mesSenNum);
            end
        end
    end
end