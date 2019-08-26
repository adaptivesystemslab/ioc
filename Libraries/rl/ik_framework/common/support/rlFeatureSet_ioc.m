classdef rlFeatureSet_ioc < rlFeatureSet % < rlModel
    properties
        % features
        time = [];
           
        q = []; % joint angle
        dq = [];
        ddq = [];
        dddq = [];
        
        x = []; % joint position
        dx = [];
        ddx = [];
        dddx = [];
        tau = []; % joint torque
        dtau = []; % torque change
        ddtau = []; % torque acceleration (Berret2011: 'effort')
        
        ek = []; % kinetic energy (modified from Berret2011's 'geodesic', which where ek = geodesic^2)
        power = []; % angular power (Berret2011: 'energy')
        
        % other params
        x_frame_name = ''; % if empty, then x/dx/ddx will be populated with the last 'body' Cart pos
        splineModel = [];
        modelStruct = [];
        
        inds_pos = [];
        inds_vel = [];
        inds_fb_pos = [];
        inds_model_pos = [];
        inds_fb_vel = [];
        inds_model_vel = [];
        
        ekfParams = [];
        
        joint_labels = [];
    end
    
    methods
        function setModelParam(obj, mdl, dt, algorithmParam, modelInstance)
            setModelParam@rlFeatureSet(obj, mdl, dt, algorithmParam);
        end
        
%         function splineModelInit(obj)
%             obj.splineModel = splineFit(enumSplineTypes.piecewiseQuintic);
%         end
%         
%         function splineModelSetViaPointsStruct(obj, splineViaPoints)
%             obj.splineModel.setViaPointsStruct(splineViaPoints);
%             obj.splineModel.calcSplineParam();
%             
%             if 0
%                 time = 0:0.01:0.6;
%                 obj.calcQFromSpline(time);
%                 
%                 figure;
%                 plot(obj.time, obj.q);
%                 hold on;
%                 for i = 1:numel(splineViaPoints)
%                     plot(splineViaPoints(i).time, splineViaPoints(i).q, 'x');
%                 end
%             end
%         end
%         
%         function splineModelSetViaPointsStructFromTrajectory(obj, t_knot)
%             [splineViaPoints.time, splineViaPoints.inds] = findClosestValue(t_knot, obj.time);
%             splineViaPoints.q = obj.q(splineViaPoints.inds, :);
%             splineViaPoints.dq = obj.dq(splineViaPoints.inds, :);
%             splineViaPoints.ddq = obj.ddq(splineViaPoints.inds, :);
%             
%             obj.splineModelSetViaPointsStruct(splineViaPoints);
%         end
%         
%         function splineModelUpdateViaPoints(obj, q_mat)
%             obj.splineModel.updateViaPoints(q_mat);
%             obj.splineModel.calcSplineParam();
%         end
% 
%         % generate q from spline data
%         function calcQFromSpline(obj, time_spline)
%             obj.time = time_spline;
%             [obj.q, obj.dq, obj.ddq] = obj.splineModel.calcQ(time_spline);
%         end
        
        % assuming EKF IK has already been done, parse the data into q
        function calcQfromEkfStates(obj, inds_fb_pos, inds_model_pos, inds_fb_vel, inds_model_vel, filter)
            if ~exist('filter', 'var')
                filter = 1;
            end
            
            if exist('inds_fb_pos', 'var')
                obj.inds_pos = [inds_fb_pos inds_model_pos];
                obj.inds_vel = [inds_fb_vel inds_model_vel];

                obj.inds_fb_pos = inds_fb_pos;
                obj.inds_model_pos = inds_model_pos;
                obj.inds_fb_vel = inds_fb_vel;
                obj.inds_model_vel = inds_model_vel;
            else
                obj.inds_pos = 1:length(obj.joint_labels);
            end
            
            if filter == 1
                q_noFilt = obj.ekf_states(:, obj.inds_pos);
                q_filt = filter_dualpassBW(q_noFilt);
                obj.q = q_filt;
            else
                obj.q = obj.ekf_states(:, obj.inds_pos);
            end

%             obj.dq = obj.ekf_states(:, obj.inds_vel);
            obj.dq = calcDerivVert(obj.q, obj.dt);
            obj.ddq = calcDerivVert(obj.dq, obj.dt);
            
%             for i = 1:length(obj.inds_pos)
%                 obj.joint_labels{i} = obj.model.joints(obj.inds_pos(i)).name;
%             end
        end
        
        % calculate all features from q
        function calcFeaturesFromQ(obj, time)
            % initialize the rest
            obj.time = time;
            obj.dddq = zeros(size(obj.q));
            
            if ~isempty(obj.algorithmParam.endEffectorName)
                obj.x = zeros(numel(obj.time), 3);
                obj.dx = zeros(numel(obj.time), 3);
                obj.ddx = zeros(numel(obj.time), 3);
                obj.dddx = zeros(numel(obj.time), 3);
            else
                obj.x = zeros(numel(obj.time), size(obj.model.bodies, 2)*3);
                obj.dx = zeros(numel(obj.time), size(obj.model.bodies, 2)*3);
                obj.ddx = zeros(numel(obj.time), size(obj.model.bodies, 2)*3);
                obj.dddx = zeros(numel(obj.time), size(obj.model.bodies, 2)*3);
            end
     
            obj.tau = zeros(size(obj.q));
            obj.dtau = zeros(size(obj.q));
            obj.ddtau = zeros(size(obj.q));
            
            obj.ek = zeros(size(obj.q));
            obj.power = zeros(size(obj.q));
            
%             zeroQ = zeros(1, size(obj.q, 2));
            
%             if obj.algorithmParam.visualize
%                 vis = rlVisualizer('vis',640,480);
%                 vis.addModel(obj.model);
%                 obj.model.forwardPosition();
%                 vis.update();
%             end

            for ind=1:numel(obj.time)
                currQ = obj.q(ind, :);
                currDQ = obj.dq(ind, :);
                currDDQ = obj.ddq(ind, :);
                
                obj.inputQDataIntoModel(currQ, currDQ, currDDQ);
                
                obj.model.inverseDynamics();
                obj.tau(ind, :) = obj.model.torque';
                
                if ~isempty(obj.algorithmParam.endEffectorName)
                    obj.x(ind, :) = ...
                        obj.model.getFrameByName(obj.algorithmParam.endEffectorName).t(1:3, 4)';
                else
                    for j = 1:numel(obj.model.bodies)
                        targetInd = (j-1)*3+1:(j-1)*3+3;
                        obj.x(ind, targetInd) = obj.model.bodies(j).t(1:3, 4)';
%                         obj.x(ind, :) = obj.model.bodies(j).t(1:3, 4)';
                    end
                end
                
                if (ind == 1 || mod(ind, 1) == 0)
                    %                     fprintf('Frame %u of %u\n', ind, numel(obj.time));
                    %
                    %                     if obj.algorithmParam.visualize
                    % %                         applyMarkersToVisualization(vis, obj.measurement_label, obj.measurement_indArray, ind);
                    %                         vis.update();
                    %                         pause(0.01);
                    %                     end
                end
                
%         f = model.getFrameByName('frame_elbow0');
%         ddx_curr = f.a(4:6) + cross(f.v(1:3),f.v(4:6)) - f.t(1:3,1:3)'*model.g;
%         ddx_curr = rotz(pi/2)*ddx_curr;
                
%         model.calculateJacobian();              model = inputData(model, q_7dof(:, i)', dq_7dof(:, i)', ddq_7dof(:, i)');
%         model.calculateJacobianDerivative();    model = inputData(model, q_7dof(:, i)', dq_7dof(:, i)', ddq_7dof(:, i)');
%         dx_curr = model.J*model.velocity;
%         dx_curr = dx_curr([1 2 3]);

%         ddx_curr1 = model.J*model.acceleration;
%         ddx_curr2 = modelr1 + ddx_curr2;
%         ddx_curr = .JdQd; %%% TODO
%         ddx_curr = ddx_curddx_curr([1 2 3]);
                
                obj.model.calculateMassMatrix();
                obj.ek(ind, :) = (obj.model.M * currDQ')'.^2;
                
%                 obj.inputQDataIntoModel(obj.q(i, :), zeroQ, zeroQ);
%                 obj.model.calculateCentrifugalCoriolis();
%                 Cdq = model.V;
%                 obj.inputQDataIntoModel(obj.q(i, :), zeroQ, zeroQ);
%                 obj.model.calculateGravity();
%                 G = model.G;
%                 obj.ep(i, :) = Cdq + G;
            end
            
            obj.dddq = calcDerivVert(obj.ddq, obj.dt);
            
            obj.dx = calcDerivVert(obj.x, obj.dt);
            obj.ddx = calcDerivVert(obj.dx, obj.dt);
            obj.dddx = calcDerivVert(obj.ddx, obj.dt);
            obj.dtau = calcDerivVert(obj.tau, obj.dt);
            obj.ddtau = calcDerivVert(obj.dtau, obj.dt);
            obj.power = obj.dq .* obj.tau;
        end
        
        function calcDdqFromTorque(obj)
            for ind=1:numel(obj.time)
                currQ = obj.q(ind, :);
                currDQ = obj.dq(ind, :);
                currDDQ = zeros(size(obj.ddq(ind, :)));
                
                obj.inputQDataIntoModel(currQ, currDQ, currDDQ);
                
                obj.model.forwardDynamics();
                obj.ddq(ind, :) = obj.model.acceleration';
            end
        end
        
        function calcTorqueFromFP(obj, fpData)
            % initialize
            obj.tau = zeros(size(obj.q));
            obj.dtau = zeros(size(obj.q));
            obj.ddtau = zeros(size(obj.q));
            
            obj.ek = zeros(size(obj.q));
            obj.power = zeros(size(obj.q));
            
            if obj.algorithmParam.visualize
                vis = rlVisualizer('vis',640,480);
                vis.addModel(obj.model);
                obj.model.forwardPosition();
                vis.update();
            end
            
            % find the two frames that corresponds with the FPs
            allFrameNames = {obj.model.bodies.name};
            indFP_R = find(ismember(allFrameNames, 'body_rankle_rballfoot') == 1);
            indFP_L = find(ismember(allFrameNames, 'body_lankle_lballfoot') == 1);
            
%             if 1
%                 obj.model.base = 'world';
%                 quse = obj.q;
%             else
                obj.model.base = 'frame_rballfoot_end';
                quse = filter_dualpassBW(obj.q);
                quse(:, 1:11) = -quse(:, 1:11);
%             end
            dquse = calcDerivVert(quse, obj.dt);
            ddquse = calcDerivVert(dquse, obj.dt);
            
            for ind=1:numel(obj.time)
%                 currQ = obj.q(ind, :);
%                 currDQ = obj.dq(ind, :);
%                 currDDQ = obj.ddq(ind, :);
                currQ = quse(ind, :);
                currDQ = dquse(ind, :);
                currDDQ = ddquse(ind, :);
                
                % update kinematics
                obj.inputQDataIntoModel(currQ, currDQ, currDDQ);
                
%                 vis.update();
%                 pause(0.01);
                
                % update external forces
%                 obj.model.bodies(indFP_R).fX = zeros(6, 1);
%                 obj.model.bodies(indFP_L).fX = zeros(6, 1);
%                 obj.model.bodies(indFP_R).fX = [fpData.M_R(ind, :) fpData.F_R(ind, :)]';
%                 obj.model.bodies(indFP_L).fX = [fpData.M_L(ind, :) fpData.F_L(ind, :)]';
                
                obj.model.inverseDynamics();
                obj.tau(ind, :) = obj.model.torque';
                
                obj.model.calculateMassMatrix();
                obj.ek(ind, :) = (obj.model.M * currDQ')'.^2;
            end
            
            obj.dtau = calcDerivVert(obj.tau, obj.dt);
            obj.ddtau = calcDerivVert(obj.dtau, obj.dt);
            obj.power = obj.dq .* obj.tau;
            
            if 1
                figure;
                t = 1:numel(obj.time);
                %             n = 4:11;
                n = 9;
                c = distinguishable_colors(length(n));
                for j = 1:length(n)
                    subplot(211);
                    toPlot = obj.q(t, n(j));
                    plot(toPlot, 'Color', c(j, :), 'DisplayName', obj.model.joints(n(j)).name);
                    
                    subplot(212);
                    toPlot = obj.tau(t, n(j));
                    plot(toPlot, 'Color', c(j, :), 'DisplayName', obj.model.joints(n(j)).name);
                    hold on
                end
                %             xlim([0 2499]);
                ylim([-1500 1500]);
                legend show
            end
        end
        
        function inputQDataIntoModel(obj, q, dq, ddq)
            obj.model.position = q;
            obj.model.velocity = dq;
            obj.model.acceleration = ddq;
            
            obj.model.forwardPosition();
            obj.model.forwardVelocity();
            obj.model.forwardAcceleration();
        end
        
        function playTrajectoryTimestep(obj, i)
            obj.model.position(:) = obj.q(i, :);
            obj.model.velocity(:) = obj.dq(i, :);
            obj.model.acceleration(:) = obj.ddq(i, :);
            
            obj.model.forwardPosition();
        end
        
        function playTrajectory(obj)
            vis = rlVisualizer('vis',640,480);
            vis.addModel(obj.model);
            
            for i=1:numel(obj.time)
                fprintf('Frame %u of %u\n', i, numel(obj.time));
                
                %Set the joint angles
                obj.playTrajectoryTimestep(i);

                vis.update();
                pause(0.01);
            end
            
            fprintf('Hit any key to continue');
            pause();
        end
        
        function loadDataFromExternal(obj, featureLoad)
            p = inputParser;
            p.KeepUnmatched = true;
            
            addOptional(p, 'time', []); % says its optional, but it's not
            addOptional(p, 'q', []);
            
            addOptional(p, 'dt', []);
            addOptional(p, 'segments', []);
            
            parse(p, featureLoad);
            
            obj.time = p.Results.time;
            obj.q = p.Results.q;
            obj.segments = p.Results.segments;
            
            if ~isempty(p.Results.dt)
                obj.dt = p.Results.dt;
            else
                obj.dt = mean(diff(obj.time));
            end
            
            obj.dq = calcDeriv(obj.q, obj.dt);
            obj.ddq = calcDeriv(obj.dq, obj.dt);
        end
        
        function saveVar = saveData(obj, filepath, baseFrameOrig, ekfTuningParam, modelInstance)
%             for i = 1:size(obj.measurement_label, 2)
%                 ind = i*3-2:i*3; % pos only
%                 mocapMarkers.(obj.measurement_label{i}) = obj.measurement_output(:, ind);
%             end
%             mocapMarkers.array = obj.measurement_output;
            
            if ~isempty(modelInstance)
                baseFrameTransform = modelInstance.model.getFrameByName(obj.model.base).t;
                
                saveVar.subjectId = modelInstance.subjectId;
                saveVar.body_height = modelInstance.body_height;
                saveVar.body_weight = modelInstance.body_weight;
                saveVar.gender = modelInstance.gender;
            else
                baseFrameTransform = eye(4);
                
                saveVar.subjectId = [];
                saveVar.body_height = [];
                saveVar.body_weight = [];
                saveVar.gender = [];
            end
            
            if ~isempty(obj.model)
                origBaseFrame = obj.model.base;
                obj.model.base = baseFrameOrig;
                obj.inputQDataIntoModel(0, 0, 0);

                for i = 1:length(obj.model.transforms)
                    modelSpecs.transforms(i).t = obj.model.transforms(i).t;
                    modelSpecs.transforms(i).name = obj.model.transforms(i).name;
                    modelSpecs.transforms(i).type = obj.model.transforms(i).type;
                    modelSpecs.transforms(i).frame_in = obj.model.transforms(i).frame_in;
                end

                for i = 1:length(obj.model.bodies)
                    modelSpecs.bodies(i).name = obj.model.bodies(i).name;
                    modelSpecs.bodies(i).mass = obj.model.bodies(i).m;
                    modelSpecs.bodies(i).com = obj.model.bodies(i).com;
                    modelSpecs.bodies(i).inertial = obj.model.bodies(i).I;
                    modelSpecs.bodies(i).length = norm(obj.model.bodies(i).t(1:3, 4));
                end
            else
                modelSpecs = [];
                origBaseFrame = '';
            end

            if ~isempty(modelInstance)
                modelStruct = modelInstance.saveModel('');
            else
                modelStruct = [];
            end

            if ~isempty(obj.x) && ~isempty(obj.model)
                for i = 1:numel(obj.model.bodies)
                    targetInd = (i-1)*3+1:(i-1)*3+3;
                    frameStr = strrep(obj.model.bodies(i).name, ':', '_');
                    endEffectorPosition.(frameStr) = obj.x(:, targetInd);
                end
            end
            endEffectorPosition.array = obj.x;
            
            for i = 1:size(obj.q, 2)
                frameStr = strrep(obj.joint_labels{i}, ':', '_');
                jointAngle.(frameStr) = obj.q(:, i);
            end
            jointAngle.array = obj.q;
            
            if ~isempty(obj.tau)
                for i = 1:size(obj.q, 2)
                    frameStr = strrep(obj.joint_labels{i}, ':', '_');
                    jointTorque.(frameStr) = obj.tau(:, i);
                end
            end
            jointTorque.array = obj.tau;
            
            saveVar.time = obj.time;
            saveVar.dt = obj.dt;
            saveVar.ekf_states = obj.ekf_states;
            
            saveVar.endEffectorPosition = endEffectorPosition;
%             saveVar.mocapMarkers = mocapMarkers;

            saveVar.jointLabels = obj.joint_labels;
            saveVar.jointAngle = jointAngle;
            saveVar.jointTorque = jointTorque;
            
            saveVar.measurementLabels = obj.measurement_labels;
            saveVar.measurementInput = obj.measurement_input;
            saveVar.measurementOutput = obj.measurement_output;
            saveVar.measurementMask = obj.measurement_mask;
            saveVar.measurementInputMatch = obj.measurement_input_match;
            saveVar.measurementOutputMatch = obj.measurement_output_match;
            saveVar.measurementEkfMatch = obj.measurement_ekf_match;

            saveVar.modelStruct = modelStruct;
            saveVar.modelSpecs = modelSpecs;
            saveVar.baseFrame = origBaseFrame;
            saveVar.baseFrameT = baseFrameTransform;
            saveVar.frameData = obj.frameData;
            saveVar.segmentData = obj.segments;
           
            saveVar.inds_pos = obj.inds_pos;
            saveVar.inds_vel = obj.inds_vel;
            saveVar.inds_fb_pos = obj.inds_fb_pos;
            saveVar.inds_model_pos = obj.inds_model_pos;
            saveVar.inds_fb_vel = obj.inds_model_vel;
            saveVar.inds_model_vel = obj.inds_model_vel;
    
            saveVar.ekf_match = obj.ekf_match;
            saveVar.ekf_matchLabels = obj.ekf_matchLabels;
            saveVar.ekfParams = ekfTuningParam;
            
            save(filepath, 'saveVar');
        end
        
        function saveVar = saveDataJointsOnly(obj, filepath, baseFrameOrig, ekfTuningParam, modelInstance)
            %             for i = 1:size(obj.measurement_label, 2)
%                 ind = i*3-2:i*3; % pos only
%                 mocapMarkers.(obj.measurement_label{i}) = obj.measurement_output(:, ind);
%             end
%             mocapMarkers.array = obj.measurement_output;

            baseFrameTransform = modelInstance.model.getFrameByName(obj.model.base).t;

%             if ~isempty(obj.model)
%                 origBaseFrame = obj.model.base;
%                 obj.model.base = baseFrameOrig;
%                 obj.inputQDataIntoModel(0, 0, 0);
% 
%                 for i = 1:length(obj.model.transforms)
%                     modelSpecs.transforms(i).t = obj.model.transforms(i).t;
%                     modelSpecs.transforms(i).name = obj.model.transforms(i).name;
%                     modelSpecs.transforms(i).type = obj.model.transforms(i).type;
%                     modelSpecs.transforms(i).frame_in = obj.model.transforms(i).frame_in;
%                 end
% 
%                 for i = 1:length(obj.model.bodies)
%                     modelSpecs.bodies(i).name = obj.model.bodies(i).name;
%                     modelSpecs.bodies(i).mass = obj.model.bodies(i).m;
%                     modelSpecs.bodies(i).com = obj.model.bodies(i).com;
%                     modelSpecs.bodies(i).inertial = obj.model.bodies(i).I;
%                     modelSpecs.bodies(i).length = norm(obj.model.bodies(i).t(1:3, 4));
%                 end
%             else
%                 modelSpecs = [];
%                 origBaseFrame = '';
%             end

            if ~isempty(modelInstance)
                modelStruct = modelInstance.saveModel('');
            else
                modelStruct = [];
            end

%             if ~isempty(obj.x) && ~isempty(obj.model)
%                 for i = 1:numel(obj.model.bodies)
%                     targetInd = (i-1)*3+1:(i-1)*3+3;
%                     endEffectorPosition.(obj.model.bodies(i).name) = obj.x(:, targetInd);
%                 end
%             end
%             endEffectorPosition.array = obj.x;
            
%             for i = 1:size(obj.q, 2)
%                 jointAngle.(obj.joint_labels{i}) = obj.q(:, i);
%             end
            jointAngle.array = obj.q;
            
%             if ~isempty(obj.tau)
%                 for i = 1:size(obj.q, 2)
%                     jointTorque.(obj.joint_labels{i}) = obj.tau(:, i);
%                 end
%             end
            jointTorque.array = obj.tau;
            
            saveVar.subjectId = modelInstance.subjectId;
            saveVar.body_height = modelInstance.body_height;
            saveVar.body_weight = modelInstance.body_weight;
            saveVar.gender = modelInstance.gender;
            
            saveVar.time = obj.time;
            saveVar.dt = obj.dt;
%             saveVar.ekf_states = obj.ekf_states;
            
%             saveVar.endEffectorPosition = endEffectorPosition;
%             saveVar.mocapMarkers = mocapMarkers;

            saveVar.jointLabels = obj.joint_labels;
            saveVar.jointAngle = jointAngle;
            saveVar.jointTorque = jointTorque;
            
%             saveVar.measurementLabels = obj.measurement_labels;
%             saveVar.measurementInput = obj.measurement_input;
%             saveVar.measurementOutput = obj.measurement_output;
%             saveVar.measurementMask = obj.measurement_mask;
%             saveVar.measurementInputMatch = obj.measurement_input_match;
%             saveVar.measurementOutputMatch = obj.measurement_output_match;
%             saveVar.measurementEkfMatch = obj.measurement_ekf_match;

            saveVar.modelStruct = modelStruct;
%             saveVar.modelSpecs = modelSpecs;
            saveVar.baseFrame = origBaseFrame;
            saveVar.baseFrameT = baseFrameTransform;
%             saveVar.frameData = obj.frameData;
%             saveVar.segmentData = obj.segments;
           
            saveVar.inds_pos = obj.inds_pos;
            saveVar.inds_vel = obj.inds_vel;
            saveVar.inds_fb_pos = obj.inds_fb_pos;
            saveVar.inds_model_pos = obj.inds_model_pos;
            saveVar.inds_fb_vel = obj.inds_model_vel;
            saveVar.inds_model_vel = obj.inds_model_vel;
    
%             saveVar.ekf_match = obj.ekf_match;
%             saveVar.ekf_matchLabels = obj.ekf_matchLabels;
            saveVar.ekfParams = ekfTuningParam;
            
            save(filepath, 'saveVar');
        end
        
        function saveVar = loadDataFromFileJointsOnly(obj, filepath, baseFrameOrig, ekfTuningParam, modelInstance)
            load(filepath);
            
            obj.time = saveVar.time;
            obj.dt = saveVar.dt;
            
%             obj.x = saveVar.endEffectorPosition.array;
            
            obj.joint_labels = saveVar.jointLabels;
            obj.q = saveVar.jointAngle.array;
            obj.tau = saveVar.jointTorque.array;
            
%             obj.measurement_labels = saveVar.measurementLabels;
%             obj.measurement_input = saveVar.measurementInput;
%             obj.measurement_output = saveVar.measurementOutput;
%             obj.measurement_mask = saveVar.measurementMask;
%             obj.measurement_ekf_match =  saveVar.measurementEkfMatch;
%             obj.measurement_input_match =  saveVar.measurementInputMatch;
%             obj.measurement_output_match =  saveVar.measurementOutputMatch;

            obj.modelStruct = saveVar.modelStruct;
%             modelSpecs = saveVar.modelSpecs;      
%             obj.baseFrame = saveVar.baseFrame;
%             obj.frameData = saveVar.frameData;
%             obj.segments = saveVar.segmentData;
            
            obj.inds_pos = saveVar.inds_pos;
            obj.inds_vel = saveVar.inds_vel;
            obj.inds_fb_pos = saveVar.inds_fb_pos;
            obj.inds_model_pos = saveVar.inds_model_pos;
            obj.inds_fb_vel = saveVar.inds_model_vel;
            obj.inds_model_vel = saveVar.inds_model_vel;
            
%             obj.ekf_match = saveVar.ekf_match;
%             obj.ekf_matchLabels = saveVar.ekf_matchLabels;
            obj.ekfParams = saveVar.ekfParams;
            
            obj.dq = calcDerivVert(obj.q, obj.dt);
            obj.ddq = calcDerivVert(obj.dq, obj.dt);
        end
        
        function saveVar = loadDataFromFile(obj, filepath)
            load(filepath);
            
            obj.time = saveVar.time;
            obj.dt = saveVar.dt;
            
            obj.x = saveVar.endEffectorPosition.array;
            
            obj.joint_labels = saveVar.jointLabels;
            obj.q = saveVar.jointAngle.array;
            obj.tau = saveVar.jointTorque.array;
            
            obj.measurement_labels = saveVar.measurementLabels;
            obj.measurement_input = saveVar.measurementInput;
            obj.measurement_output = saveVar.measurementOutput;
            obj.measurement_mask = saveVar.measurementMask;
            obj.measurement_ekf_match =  saveVar.measurementEkfMatch;
            obj.measurement_input_match =  saveVar.measurementInputMatch;
            obj.measurement_output_match =  saveVar.measurementOutputMatch;

            obj.modelStruct = saveVar.modelStruct;
%             modelSpecs = saveVar.modelSpecs;      
            obj.baseFrame = saveVar.baseFrame;
            obj.frameData = saveVar.frameData;
            obj.segments = saveVar.segmentData;
            
            obj.inds_pos = saveVar.inds_pos;
            obj.inds_vel = saveVar.inds_vel;
            obj.inds_fb_pos = saveVar.inds_fb_pos;
            obj.inds_model_pos = saveVar.inds_model_pos;
            obj.inds_fb_vel = saveVar.inds_model_vel;
            obj.inds_model_vel = saveVar.inds_model_vel;
            
            obj.ekf_match = saveVar.ekf_match;
            obj.ekf_matchLabels = saveVar.ekf_matchLabels;
            obj.ekfParams = saveVar.ekfParams;
            
            obj.dq = calcDerivVert(obj.q, obj.dt);
            obj.ddq = calcDerivVert(obj.dq, obj.dt);
        end
        
        function cropLoadedData(obj, cropInds)
            obj.time = obj.time(cropInds);
            obj.q = obj.q(cropInds, :);
%             obj.x = obj.x(cropInds, :);
%             obj.tau = obj.tau(cropInds, :);
            
            if ~isempty(obj.segments)
                ind = 0;
                for i = 1:length(obj.segments)
                    if obj.segments(i).timeStart > obj.time(1) && obj.segments(i).timeEnd > obj.time(1) && ...
                            obj.segments(i).timeStart < obj.time(end) && obj.segments(i).timeEnd < obj.time(end)
                        ind = ind + 1;
                        segmentsKeep(ind) = obj.segments(i);
                    end
                end
                obj.segments = segmentsKeep;
            end
        end
        
        function figureMatchMatrix(obj, filepath, modelInstance)
%             indsToLoad{1} = modelInstance.inds_mocapMarkers;
            indsToLoad{1} = 1:size(obj.ekf_match, 2);
            colours = distinguishable_colors(size(obj.ekf_match, 2));
            
            for j = 1:length(indsToLoad)
                currInds = indsToLoad{j};
                h = figure('position', 1.0e+03 *[ 0.0010    0.0410    1.5360    0.7488]);
                hold on;    
                for i = currInds
                    ekfInd = obj.ekf_match(1, i);
                    if ekfInd > 0
                        plot(obj.time, obj.ekf_match(:, i), 'Color', colours(i, :), 'DisplayName', obj.ekf_matchLabels{ekfInd});
                    end
                end
                
                yticks(1:length(obj.ekf_matchLabels));
                yticklabels(obj.ekf_matchLabels);
                xlabel('Time [s]');
                ylabel('Input Mes Label');
                
                if exist('filepath', 'var')
                    saveas(h, [filepath '_' num2str(j)], 'png');
                    saveas(h, [filepath '_' num2str(j)], 'fig');
                end
                close(h);
            end
        end
        
        function figureFrameData(obj, filepath, modelInstance)
            switch length(obj.frameData)
                case 1
                    
                case 2
                    
                case 3
                    
            end
        end
        
        function figureJointData(obj, filepath, modelInstance, plotType)
            if ~exist('plotType', 'var')
                plotType = 'angle'; % velocity angleInd
            end
            
            lenJointLabels = length(obj.joint_labels);
            
            if isempty(obj.time)
                obj.time = obj.dt*(0:size(obj.q, 1)-1);
            end
            
            % figure out which joints go together
            switch plotType
                case {'angle', 'velocity'}
                    currNameSplit = strsplit(obj.joint_labels{1}(1:end), '_');
                    currName = currNameSplit{2};
                    
                    nameArray = {};
                    nameArray{1, 1} = currName;
                    nameArray{1, 2} = 1;
                    
                    for i = 2:lenJointLabels
                        currNameSplit = strsplit(obj.joint_labels{i}(1:end), '_');
                        currName = currNameSplit{2};
                        
                        % normal case
                        findInd = find(strcmpi(currName, nameArray(:, 1)));
                        if ~isempty(findInd)
                            nameArray{findInd, 2} = [nameArray{findInd, 2} i];
                        else
                            nameArray{size(nameArray, 1)+1, 1} = currName;
                            nameArray{end, 2} = i;
                        end
                    end
                    
                case 'angleInd'
                    nameArray = {};
                    
                    for i = 1:lenJointLabels
                        currNameSplit = strsplit(obj.joint_labels{1}(1:end), '_');
                        currName = currNameSplit{2};
                        
                        nameArray{i, 1} = currName;
                        nameArray{i, 2} = i;
                    end
            end
          
            h = figure('position', 1.0e+03 *[ 0.0010    0.0410    1.5360    0.7488]);
            associatedJoints = size(nameArray, 1);
            for i = 1:associatedJoints
                switch associatedJoints
                    case {1, 2, 3, 4}
                        ax(i) = subplot(2, 2, i);
                    case {5, 6}
                        ax(i) = subplot(2, 3, i);
                    case {7, 8, 9, 11, 12}
                        ax(i) = subplot(3, 4, i);
                    case {13, 14, 15}
                        ax(i) = subplot(3, 5, i);
                    case {16, 17, 18, 19, 20, 21}
                        ax(i) = subplot(3, 7, i);
                    otherwise
                        ax(i) = subplot(4, 9, i);
                end
                
                currName = nameArray{i, 1};
                currJoints = nameArray{i, 2};
                currJointLabels = length(nameArray{i, 2});
                
                colours = distinguishable_colors(currJointLabels);
                
                hold on;
                for j = 1:length(currJoints)
                    currJointName = obj.joint_labels{currJoints(j)};
                    
                    if ~isempty(modelInstance.jointRangeOfMotion)  
                        for k = 1:length(modelInstance.jointRangeOfMotion)
                            if strcmpi(modelInstance.jointRangeOfMotion{k}{1}, currJointName)
                                lb = (modelInstance.jointRangeOfMotion{k}{2}(1)-modelInstance.jointRangeOfMotionTolerance);
                                ub = (modelInstance.jointRangeOfMotion{k}{2}(2)+modelInstance.jointRangeOfMotionTolerance);
                                                                
                                tLen = length(obj.time);
                                ydata = [lb ub];
                                
%                                 xdata = [obj.time(floor(tLen/3)) obj.time(floor(2*tLen/3))];
%                                 rectangle('Position', [-1 ydata(1) xdata(2)*2 ydata(2)-ydata(1)], ...
%                                     'FaceColor', [153,216,201]/255, 'LineStyle', 'none');
                                
                                xdata = [obj.time(1) obj.time(end)];
                                ydata = ones(2, 1)*lb;
                                plot(xdata, ydata, 'Color', colours(j, :), 'LineStyle', '--',  ...
                                    'DisplayName', [currJointName '_lim']);
                                
                                ydata = ones(2, 1)*ub;
                                plot(xdata, ydata, 'Color', colours(j, :), 'LineStyle', '--',  ...
                                    'DisplayName', [currJointName '_lim']);
                            end
                        end
                    end

                    switch plotType
                        case 'angle'
                            plot(obj.time, rad2deg(obj.q(:, currJoints(j))), 'Color', colours(j, :), ...
                                'DisplayName', currJointName);
                            ylabel('[deg]');
                            
                        case 'velocity'
                            plot(obj.time, rad2deg(obj.dq(:, currJoints(j))), 'Color', colours(j, :), ...
                                'DisplayName', currJointName);
                            ylabel('[deg/s]');
                    end
                end
                
                title(currName);
                
                if ~isempty(obj.segments)
                    %                     [minY, maxY, h_line] = plotBoxes(segmentInfo, colorToUse, offset, minY, maxY)
                    plotBoxes(obj.segments, [1 0 0], 0, -7*pi/8, 7*pi/8);
                end
            end
            
            linkaxes(ax, 'x');
            xlim([obj.time(1) obj.time(end)])
            
%             legend show;
            
            if exist('filepath', 'var')
                saveas(h, filepath, 'png');
                saveas(h, filepath, 'fig');
            end
            close(h);
        end
        
        function figureMeasurements(obj, filepath, ekf_markerMask, ekf_eventhandler, ekf_markerTally)
            plot_measurements(obj, obj, obj.measurement_input, obj.measurement_output, ...
                'obs', 'est', filepath, ekf_markerMask, ekf_eventhandler, ekf_markerTally);
        end
        
        function resample(obj, newDt)
            if newDt == obj.dt
                return;
            end
            
            newTime = obj.time(1):newDt:obj.time(end);
            qNew = interp1(obj.time, obj.q, newTime);
            
            obj.q = qNew;
            obj.dq = calcDerivVert(obj.q, obj.dt);
            obj.ddq = calcDerivVert(obj.dq, obj.dt);
            
            obj.dt = newDt;
        end
    end
end