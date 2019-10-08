classdef IOCInstance < handle
    properties
        delta = 5e-4; % delta for numerical gradients
        winSize = 3; % window length to consider to calculate the basis features
        
        dynamicModel; 
        features; % cost function features
        basisFeatureFlags; % a list of the basis features to calculate 
        dt;
        
        jointRequiredStr;
        jointRequiredInd;
        frameRequiredStr;
    end
    
    methods
        function obj = IOCInstance(dynamicModel, dt)
            % constructor
            obj.dynamicModel = dynamicModel;
            obj.dt = dt;
        end
        
        function init(obj, jsonBlob)
            % translate the json file content into feature information to
            % be extracted from the dynamic model
            rlModel = obj.dynamicModel.model;
            
            % initialize the variables to something temp
            obj.features = featuresCalc();
            basisFeatures = basisFeaturesEnums.dddq;
            jointRequired = {};
            frameRequired = {};
            
            % initialize the feature list
            for i = 1:length(jsonBlob.candidateFeatures)
                currJsonBlob = jsonBlob.candidateFeatures(i);
                obj.features(i) = featuresCalc(currJsonBlob);
                
                basisFeatures = [basisFeatures obj.features(i).bf];

                switch obj.features(i).jointType
                    case 'joint' 
                        jointRequired = [jointRequired obj.features(i).param];
                        
                    case 'frame'
                        frameRequired = [frameRequired obj.features(i).param];
                end
            end
            
            % update the full basis feature list so calcFeatures will know
            % which basis features are required
            basisFeatures = unique(basisFeatures(2:end)); % remove the initial one
            allBasisFeatures = enumeration(basisFeatures(1));
            for i = 1:length(allBasisFeatures)
                featureUse = find(ismember(allBasisFeatures(i), basisFeatures));
                if isempty(featureUse)
                    featureUse = 0;
                end
                obj.basisFeatureFlags.(char(allBasisFeatures(i))) = featureUse;
            end
            
            % make a subindex for joint features wrt to the state length,
            % and update the individual cost function parameters with it
            obj.jointRequiredStr = unique(jointRequired);
            obj.jointRequiredStr = obj.jointRequiredStr(:); % ensure dim size
            obj.jointRequiredInd = zeros(size(obj.jointRequiredStr, 2), 1);
            for i = 1:length(obj.jointRequiredStr)
                obj.jointRequiredInd(i) = find(ismember({rlModel.joints.name}, obj.jointRequiredStr{i}));
            end
            
            % make a subindex for frame features
            obj.frameRequiredStr = unique(frameRequired);
            obj.frameRequiredStr = obj.frameRequiredStr(:); % ensure dim size
            
            % now assign these parameters back into the feature listings
            for i = 1:length(obj.features)     
                switch obj.features(i).jointType
                    case 'joint'
                        for j = 1:length(obj.features(i).param)
                            obj.features(i).subIndices(j) = find(ismember(obj.jointRequiredStr, obj.features(i).param(j)));
                        end
                        
                    case 'frame'
                        for j = 1:length(obj.features(i).param)
                            inds = (1:3)+(j-1)*3;
                            obj.features(i).subIndices(inds) = find(ismember(obj.frameRequiredStr, obj.features(i).param(j)));
                        end
                end
            end
            
            obj.dynamicModel.addEndEffectors(obj.frameRequiredStr);
        end
        
        function featureVals = calcFeatures(obj, state, control)
            % calculate the required basis features, then determine the cost 
            % function features determined by obj.features. the features 
            % are returned as a matrix with time as rows and features as
            % columns
            lenTime = size(state, 1);   
            lenFeatures = length(obj.features);
            lenJointRequired = size(obj.jointRequiredStr, 1);
            lenFrameRequired = size(obj.frameRequiredStr, 1)*3;
            
            % decoding the input, in the order as passed in
            [qState, dqState] = decodeState(state);
            tauControl = control;
            
            % initializing the basis features variables and rearranging
            % dq/tau if needed, in the jointRequired order
            if obj.basisFeatureFlags.dq
                dq = dqState(:, obj.jointRequiredInd);
            end
            
            if obj.basisFeatureFlags.tau
                tau = tauControl(:, obj.jointRequiredInd);
            end
            
            if obj.basisFeatureFlags.ddq
                ddq = zeros(lenTime, lenJointRequired);
            end
            
            if obj.basisFeatureFlags.kinetic
                kinetic = zeros(lenTime, lenJointRequired);
            end
            
            if obj.basisFeatureFlags.ddx
                ddx = zeros(lenTime, lenFrameRequired);
            end

            for i = 1:lenTime
                currQ = qState(i, :);
                currDq = dqState(i, :);
                currTau = tauControl(i, :);
                obj.dynamicModel.updateState(currQ, currDq);
                
                if obj.basisFeatureFlags.ddq
                    ddqTemp = obj.dynamicModel.forwardDynamics(currTau, 1);
                    ddq(i, :) = ddqTemp(:, obj.jointRequiredInd);
                end
                
                if obj.basisFeatureFlags.kinetic
                    obj.dynamicModel.model.calculateMassMatrix();
                    kineticTemp = (currDq * obj.dynamicModel.model.M * currDq')';
                    kinetic(i, :) = kineticTemp(:, obj.jointRequiredInd);
                end
                
                if obj.basisFeatureFlags.ddx
                    for j = 1:length(obj.frameRequiredStr)
                        inds = (1:3)+(j-1)*3;
                        ddxTemp =  obj.dynamicModel.getEndEffectorAcceleration(ddqTemp, j);
                        ddx(i, inds) = ddxTemp;
                    end
                end
            end
            
            if obj.basisFeatureFlags.dddq
                if lenTime > 1
                    dddq = calcDerivVert(ddq, obj.dt);
                else
                    dddq = zeros(lenTime, lenJointRequired);
                end
            end
            
            if obj.basisFeatureFlags.dddx
                if lenTime > 1
                    dddx = calcDerivVert(ddx, obj.dt);
                else
                    dddx = zeros(lenTime, lenFrameRequired);
                end
            end
            
            if obj.basisFeatureFlags.dtau
                if lenTime > 1
                    dtau = calcDerivVert(tau, obj.dt);
                else
                    dtau = zeros(lenTime, lenJointRequired);
                end
            end
            
            if obj.basisFeatureFlags.ddtau
                if lenTime > 2
                    ddtau = calcDerivVert(dtau, obj.dt);
                else
                    ddtau = zeros(lenTime, lenJointRequired);
                end
            end
            
            % now calculate the cost function features
            featureVals = zeros(lenTime, lenFeatures);
            
            % assemble the basis features at timestep i, since IOC only
            % considers the current timestep
            if obj.basisFeatureFlags.dq
                basisFeatureVals.dq = dq;
            end
            
            if obj.basisFeatureFlags.ddq
                basisFeatureVals.ddq = ddq;
            end
            
            if obj.basisFeatureFlags.dddq
                basisFeatureVals.dddq = dddq;
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
            
            if obj.basisFeatureFlags.kinetic
                basisFeatureVals.kinetic = kinetic;
            end
            
            for j = 1:lenFeatures
                % calculate the cost function features, where the i is
                % time and j is the feature index
                featureVals(:, j) = obj.features(j).calcFeature(basisFeatureVals);
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
                obj.dynamicModel.updateState(currQ, currDq);
                
                ddq(i, :) = obj.dynamicModel.forwardDynamics(currTau, 1);       
            end
            
            % calculating q_(t+1)
            qtp1 = q + dq*obj.dt + ddq*(obj.dt^2)/2;
            dqtp1 = dq + ddq*obj.dt;
            
%             dynVals = [qtp1 dqtp1];
            dynVals = zeros(lenTime, lenState);
            dynVals(:, 1:2:lenState) = qtp1;
            dynVals(:, 2:2:lenState) = dqtp1;
        end
        
        function gradient = getDynamicsStateDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);            
           
            gradient = zeros(size(stateOrig, 2), size(stateOrig, 2));
            for i = lenTime % preturb only the last timestep
                for j = 1:size(stateOrig, 2) % iterate through the states        
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    state(end, j) = state(end, j) + obj.delta; % add perturb to a single state
                    cf1full = obj.calcDynamics(state, control); % calculate the dynamics after perturbing
                    cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in
                    
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    state(end, j) = state(end, j) - obj.delta;
                    cf2full = obj.calcDynamics(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf1 - cf2) / (2*obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function gradient = getDynamicsControlDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);
            
            gradient = zeros(size(stateOrig, 2), size(controlOrig, 2));
            for i = lenTime % preturb only the last timestep
                for j = 1:size(controlOrig, 2) % iterate through the pos states
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    control(end, j) = control(end, j) + obj.delta; % add perturb to a single control
                    cf1full = obj.calcDynamics(state, control); % calculate the dynamics after perturbing
                    cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in
                    
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    control(end, j) = control(end, j) - obj.delta;
                    cf2full = obj.calcDynamics(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf1 - cf2) / (2*obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function gradient = getFeaturesStateDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);
            
            gradient = zeros(size(obj.features, 2), size(stateOrig, 2));
            for i = lenTime % preturb only the last timestep
                for j = 1:size(stateOrig, 2) % iterate through the pos states
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    state(end, j) = state(end, j) + obj.delta; % add perturb to a single state
                    cf1full = obj.calcFeatures(state, control); % calculate the features after perturbing
                    cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in
                    
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    state(end, j) = state(end, j) - obj.delta;
                    cf2full = obj.calcFeatures(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf1 - cf2) / (2*obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function gradient = getFeaturesControlDerivatives(obj, stateOrig, controlOrig)
            % use central difference to calculate the gradient
            lenTime = size(stateOrig, 1);
            
            gradient = zeros(size(obj.features, 2), size(controlOrig, 2));
            for i = lenTime % preturb only the last timestep
                for j = 1:size(controlOrig, 2) % iterate through the pos states
                    % crop the data being considered by
                    % calcFeatures/calcDynamics so it doesn't need to
                    % iterate through the whole time array
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    control(end, j) = control(end, j) + obj.delta; % add perturb to a single control
                    cf1full = obj.calcFeatures(state, control); % calculate the features after perturbing
                    cf1 = cf1full(end, :); % keep only timestep i of the returned cf, not the rest of the window passed in
                    
                    [state, control] = setDataLength(stateOrig, controlOrig, obj.winSize, i);
                    control(end, j) = control(end, j) - obj.delta;
                    cf2full = obj.calcFeatures(state, control);
                    cf2 = cf2full(end, :);
                    
                    gradient(:, j) = (cf1 - cf2) / (2*obj.delta); % calculate the numerical gradient
                end
            end
        end
        
        function [df_dx, df_du, dp_dx, dp_du] = getDerivativesNewObservation(obj, control, state)
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
            dp_du = dp_du';
        end
    end
end