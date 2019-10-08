classdef IOCInstanceNew < handle
    properties(Constant)
        delta = 1e-6; % delta for numerical gradients
        winSize = 3; % window length to consider to calculate the basis features
        
        filterFeature = 0; % if enabled, winsize should be at least 10, to apply a weak LPF to reduce high freq jitters
    end
    
    properties
        filtFreq = 0.25;
        filtSamp = 0;
        filtOrder = 1;
        
        dynamicModel; 
        features; % cost function features
        basisFeatureFlags; % a list of the basis features to calculate 
        dt;
        
        fullJointNames;
        fullFrameNames;
    end
    
    methods
        function obj = IOCInstanceNew(dynamicModel, dt)
            % constructor
            obj.dynamicModel = dynamicModel;
            obj.dt = dt;
        end
        
        function init(obj, jsonBlob)
            % translate the json file content into feature information to
            % be extracted from the dynamic model
            obj.dynamicModel.forwardKinematics();
            rlModel = obj.dynamicModel.model;
            
            % initialize the variables to something temp
            obj.features = featuresCalc();
            basisFeatures = basisFeaturesEnums.q;
            frameRequired = {};
            
            % init the json file and create the initial list of all joints
            %  and frames that will be considered
            for i = 1:length(jsonBlob.candidateFeatures)
                currJsonBlob = jsonBlob.candidateFeatures(i);
                obj.features(i) = featuresCalc(currJsonBlob);
                frameRequired = [frameRequired; obj.features(i).frameNames];
                
                % specific cases that requires frame information to be
                % added
                switch obj.features(i).feature
                    case featuresEnums.cartDisplacementSumSqu
                        frameRequired = [frameRequired; obj.features(i).refFrameNames];
                end
            end
            
            % use the model state order for the joints for consistency
            obj.fullJointNames = {rlModel.joints.name};
            obj.fullJointNames = obj.fullJointNames(:);
            
            % use the alphabetical frame order for joint pos
            obj.fullFrameNames = unique(frameRequired);
            obj.dynamicModel.addEndEffectors(obj.fullFrameNames);
            obj.dynamicModel.forwardKinematics();
            
            % now assign these parameters back into the feature listings
            for i = 1:length(obj.features)     
                obj.features(i).initialize(obj.fullJointNames, obj.fullFrameNames, obj.dynamicModel);
                basisFeatures = [basisFeatures obj.features(i).bf];
            end
            
            % update the full basis feature list so calcFeatures will know
            % which basis features are required
            basisFeatures = unique(basisFeatures(2:end)); % remove the initial one, q
            allBasisFeatures = enumeration(basisFeatures(1));
            for i = 1:length(allBasisFeatures)
                featureUse = find(ismember(allBasisFeatures(i), basisFeatures));
                
                if isempty(featureUse)
                    featureUse = 0;
                end
                
                obj.basisFeatureFlags.(char(allBasisFeatures(i))) = featureUse;
            end
        end
        
        function featureVals = calcFeatures(obj, state, control)
            % calculate the required basis features, then determine the cost 
            % function features determined by obj.features. the features 
            % are returned as a matrix with time as rows and features as
            % columns
            lenTime = size(state, 1);
            lenFeatures = length(obj.features);
            lenFullJoint = size(obj.fullJointNames, 1);
            lenFullFrame = size(obj.fullFrameNames, 1)*3;
            
            % decoding the input, in the order as passed in
            [q, dq] = decodeState(state);
            tau = control;
            
            % initializing the basis features variables 
            if obj.basisFeatureFlags.ddq
                ddq = zeros(lenTime, lenFullJoint);
            end

            if obj.basisFeatureFlags.x
                x = zeros(lenTime, lenFullFrame);
            end
            
            if obj.basisFeatureFlags.dx
                dx = zeros(lenTime, lenFullFrame);
            end
            
            if obj.basisFeatureFlags.ddx
                ddx = zeros(lenTime, lenFullFrame);
            end
            
            if obj.basisFeatureFlags.M
                M = zeros(lenFullJoint, lenFullJoint, lenTime);
            end
            
            if obj.basisFeatureFlags.R
                R = zeros(3, lenFullFrame * 3, lenTime); % 3 rows x 3*numFrames x timesteps
            end

            for i = 1:lenTime
                currQ = q(i, :);
                currDq = dq(i, :);
                currTau = tau(i, :);
                obj.dynamicModel.updateState(currQ, currDq);
                
                if obj.basisFeatureFlags.ddq
                    ddq(i, :) = obj.dynamicModel.forwardDynamicsQDqTau(currQ, currDq, currTau);
                end
                
                if obj.basisFeatureFlags.x
                    for j = 1:length(obj.fullFrameNames)
                        inds = (1:3)+(j-1)*3;
                        x(i, inds) = obj.dynamicModel.getEndEffectorPosition(j);
                    end
                end
                
                if obj.basisFeatureFlags.R
                    for j = 1:length(obj.fullFrameNames)
                        inds = (1:3)+(j-1)*3;
                        R(:, inds, i) = obj.dynamicModel.getEndEffectorPosition(j);
                    end
                end
                
                if obj.basisFeatureFlags.dx
                    for j = 1:length(obj.fullFrameNames)
                        inds = (1:3)+(j-1)*3;
                        dx(i, inds) = obj.dynamicModel.getEndEffectorVelocity(j);
                    end
                end
                
                if obj.basisFeatureFlags.ddx
                    for j = 1:length(obj.fullFrameNames)
                        inds = (1:3)+(j-1)*3;
                        ddx(i, inds) = obj.dynamicModel.getEndEffectorAcceleration(ddq(i, :), j);
                    end
                end
                
                if obj.basisFeatureFlags.M
                    obj.dynamicModel.updateState(currQ, currDq);
                    obj.dynamicModel.model.calculateMassMatrix();
                    M(:, :, i) = obj.dynamicModel.model.M; 
                end
            end
            
            % for features that need to be numerically diff
            if obj.basisFeatureFlags.dddq
                minWinLen = 1;
                lenWidth = lenFullJoint;
                dddq = obj.numDiff(ddq, minWinLen, lenTime, lenWidth);
            end
            
            if obj.basisFeatureFlags.dddx
                minWinLen = 1;
                lenWidth = lenFullFrame;
                dddx = obj.numDiff(ddx, minWinLen, lenTime, lenWidth);
            end
            
            if obj.basisFeatureFlags.dtau
                minWinLen = 1;
                lenWidth = lenFullJoint;
                dtau = obj.numDiff(tau, minWinLen, lenTime, lenWidth);
            end
            
            if obj.basisFeatureFlags.ddtau
                minWinLen = 2;
                lenWidth = lenFullJoint;
                ddtau = obj.numDiff(dtau, minWinLen, lenTime, lenWidth);
            end
            
            % now calculate the cost function features
            featureVals = zeros(lenTime, lenFeatures);
            
            % assemble the basis features
            if obj.basisFeatureFlags.q
                basisFeatureVals.q = q;
            end
            
            if obj.basisFeatureFlags.dq
                basisFeatureVals.dq = dq;
            end
            
            if obj.basisFeatureFlags.ddq
                basisFeatureVals.ddq = ddq;
            end
            
            if obj.basisFeatureFlags.dddq
                basisFeatureVals.dddq = dddq;
            end
            
            if obj.basisFeatureFlags.x
                basisFeatureVals.x = x;
            end
            
            if obj.basisFeatureFlags.dx
                basisFeatureVals.dx = dx;
            end
            
            if obj.basisFeatureFlags.ddx
                basisFeatureVals.ddx = ddx;
            end
            
            if obj.basisFeatureFlags.dddx
                basisFeatureVals.dddx = dddx;
            end
            
            if obj.basisFeatureFlags.tau
                basisFeatureVals.tau = tau;
            end
            
            if obj.basisFeatureFlags.dtau
                basisFeatureVals.dtau = dtau;
            end
            
            if obj.basisFeatureFlags.ddtau
                basisFeatureVals.ddtau = ddtau;
            end
            
            if obj.basisFeatureFlags.M
                basisFeatureVals.M = M;
            end

            for j = 1:lenFeatures
                % calculate the cost function features, where the i is
                % time and j is the feature index
                featureVals(:, j) = obj.features(j).calcFeature(basisFeatureVals);
            end
        end
        
        function filtData = numDiff(obj, preDiff, minWinLen, lenTime, lenWidth)
            if obj.filterFeature && lenTime >= obj.winSize
                % want to filter + window large enough to filter
                rawData = calcDerivVert(preDiff, obj.dt);
                filtDataTemp = filter_bw_lpf(rawData(2:end, :), obj.filtFreq, obj.filtOrder); % don't include the zero in the front from the diff
                filtData = [zeros(1, lenWidth); ... 
                    filtDataTemp]; % insert the zero back in so the vector is the correct width
            elseif lenTime >= minWinLen
                % window not large enough to filter
                filtData = calcDerivVert(preDiff, obj.dt);
            else
                % window not large enough to calc this feature
                filtData = zeros(lenTime, lenWidth);
            end
            
            if 0
                
                clf
                if size(rawData, 2) == 3
                    for j = 1:3
                        subplot(1, 3, j);
                        plot(1:lenTime, rawData(:, j), 'o'); hold on;
                        plot(2:lenTime, filtData(2:end, j));
                    end
                else
                    for j = 1:15
                        subplot(3, 5, j);
                        plot(1:lenTime, rawData(:, j), 'o'); hold on;
                        plot(2:lenTime, filtData(2:end, j));
                    end
                end
                
            end
        end
        
        function setFeatureNormalization(obj, state, control)
            lenFeatures = length(obj.features);
            
            f = obj.calcFeatures(state, control);
            for i = 1:lenFeatures
                if obj.features(i).normCoeff == 1
                    % if default normalization, set it
                    coeff = mean(f(:, i));
                    obj.features(i).normCoeff = coeff;
                end
            end
        end
        
        function dynVals = calcDynamics(obj, state, control)
%             lenDofs = size(state, 2)/2;
            lenState = size(state, 2);
            lenTime = size(state, 1);
            
            [q, dq] = decodeState(state);
            tau = control;
            
            ddq = zeros(size(q));
            
            % ddq
            for i = 1:lenTime
                currQ = q(i, :);
                currDq = dq(i, :);
                currTau = tau(i, :);
                
                ddq(i, :) = obj.dynamicModel.forwardDynamicsQDqTau(currQ, currDq, currTau);
                
%                 obj.dynamicModel.updateState(currQ, currDq);
%                 ddq(i, :) = obj.dynamicModel.forwardDynamics(currTau, 1);       
            end
            
            % calculating q_(t+1)
%             qtp1   = q +  dq*obj.dt + ddq*(0.5*obj.dt^2);
            qtp1   = q +  dq*obj.dt;
            dqtp1 = dq + ddq*obj.dt;
            
            dynVals = encodeState(qtp1, dqtp1);
        end
        
        function gradient = getDynamicsStateDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);            
           
            gradient = zeros(size(stateOrig, 2), size(stateOrig, 2));
            for i = lenTime % preturb only the last timestep
                % don't need to perturb things
                [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                cf1full = obj.calcDynamics(state, control); % calculate the dynamics after perturbing
                cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in
                
                for j = 1:size(stateOrig, 2) % iterate through the states        
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    state(end, j) = state(end, j) + obj.delta;
                    cf2full = obj.calcDynamics(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf2 - cf1) / (obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function gradient = getDynamicsControlDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);
            
            gradient = zeros(size(stateOrig, 2), size(controlOrig, 2));
            for i = lenTime % preturb only the last timestep
                [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                cf1full = obj.calcDynamics(state, control); % calculate the dynamics after perturbing
                cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in

                for j = 1:size(controlOrig, 2) % iterate through the pos states
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    control(end, j) = control(end, j) + obj.delta;
                    cf2full = obj.calcDynamics(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf2 - cf1) / (obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function gradient = getFeaturesStateDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);
            
            gradient = zeros(size(obj.features, 2), size(stateOrig, 2));
            for i = lenTime % preturb only the last timestep
                [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                cf1full = obj.calcFeatures(state, control); % calculate the features after perturbing
                cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in

                for j = 1:size(stateOrig, 2) % iterate through the pos states
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array                    
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    state(end, j) = state(end, j) + obj.delta;
                    cf2full = obj.calcFeatures(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf2 - cf1) / (obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function gradient = getFeaturesControlDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);
            
            gradient = zeros(size(obj.features, 2), size(controlOrig, 2));
            for i = lenTime % preturb only the last timestep
                [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                cf1full = obj.calcFeatures(state, control); % calculate the features after perturbing
                cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in
                
                for j = 1:size(controlOrig, 2) % iterate through the pos states
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    control(end, j) = control(end, j) + obj.delta;
                    cf2full = obj.calcFeatures(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf2 - cf1) / (obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function [df_dx, df_du, dp_dx, dp_du] = getDerivativesNewObservation(obj, state, control)
            % Get derivatives of dynamics wrt to state vector
            df_dx = obj.getDynamicsStateDerivatives(state, control);

            % Get derivatives of dynamics wrt to control signal
            df_du = obj.getDynamicsControlDerivatives(state, control);
            
            % Get derivatives of cost functions wrt to state vector
            dp_dx = obj.getFeaturesStateDerivatives(state, control);

            % Get derivatives of cost functions wrt to control signal
            dp_du = obj.getFeaturesControlDerivatives(state, control);
            
            % apply transpose and state rearrangement to meet WJ's notation
%             df_dx = df_dx';
%             df_du = df_du';
%             dp_dx = dp_dx';
%             dp_du = dp_du';
        end
        
        function featureList = getFeatureListAsString(obj)
            for i = 1:length(obj.features)
                featureList{i} = char(obj.features(i).feature);
            end
        end
    end
end