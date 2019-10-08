classdef featuresCalc < handle
    % denotes the type of features to calculate from the dynamic model.
    % this class is used in conjunction with IOCInstance
    properties
        % loaded feature type from .json. this is a enum of featureEnums
        feature
        % loaded parameter form .json. these are modifiers that are parsed 
        % by this function and IOCInstance to determine how to calculate 
        % the cost function
        jointNames % which joint frame to use
        frameNames % which transform frame to use
        bodyNames % if calculating COM and want true COM, the corresponding body names to the frames needs to be loaded if weights are not manually defined
        refFrameNames % reference frame or pose
        weights % weights for cfs that require it, or determine xyz for bounding box
        
        % determines if the frameInds belongs to frames (ie cartesian) or
        % joints (ie revolute)
        jointType
        % frameInds of the full array used. this is calculated by
        % IOCInstance from obj.param
        jointInds
        frameInds
        referenceInds
        referenceVals
        referenceTimes
        windowLength
        
        normCoeff
        
        % denotes all the basis features that the cost feature requires,
        % including any time dependencies on state/control variables
        bf 
        % denotes the cost feature anon function used to calculate this
        % feature
        cf
    end
    
    methods
        function obj = featuresCalc(jsonBlob)
            % to add a new cost function, add an entry in featuresEnums
            % first, then add a corresponding entry in the switch statement
            % below
            
            if nargin == 0 % initializing an instance from the outside
                return;
            end

            parseJsonFeatureList(obj, jsonBlob);
        end
        
        function parseJsonFeatureList(obj, jsonBlob)
            if iscell(jsonBlob)
                % if the features in the json file do not have the same
                % field names, it may be created as a cell array, so we
                % need to break it out of the cell array before use
                jsonBlob = jsonBlob{1};
            end
            
            % uniform parsing
            externalParamParser = inputParser;
            externalParamParser.KeepUnmatched = true;
            addOptional(externalParamParser, 'feature', []);
            addOptional(externalParamParser, 'jointNames', []);
            addOptional(externalParamParser, 'frameNames', []);
            addOptional(externalParamParser, 'bodyNames', []);
            addOptional(externalParamParser, 'weights', []);
            addOptional(externalParamParser, 'refFrameNames', []);            
            addOptional(externalParamParser, 'windowLength', []);
            addOptional(externalParamParser, 'normCoeff', 1);
            addOptional(externalParamParser, 'refTimes', []);

            parse(externalParamParser, jsonBlob);
            jsonParsed = externalParamParser.Results;
            
            obj.feature = featuresEnums.(jsonParsed.feature); % enum parse
            obj.jointNames = jsonParsed.jointNames(:);
            obj.frameNames = jsonParsed.frameNames(:);
            obj.bodyNames = jsonParsed.bodyNames(:);
            obj.refFrameNames = jsonParsed.refFrameNames;
            obj.referenceTimes = jsonParsed.refTimes(:);
            
            obj.weights = jsonParsed.weights;
            obj.windowLength = jsonParsed.windowLength;
            
            obj.normCoeff = jsonParsed.normCoeff;
        end
        
        function initialize(obj, fullJointNames, fullFrameNames, dynamicModel)
            sumSqu = @(x) sum(x.^2, 2);
%             sumAbs = @(x) sum(abs(x), 2);
%             maxAbs = @(x) max(abs(x), [], 2);
            
            obj.bf = basisFeaturesEnums.dddx;  % initializing to something temp for now
            
            obj.jointInds = featuresCalc.setJointInds(fullJointNames, obj.jointNames);
            obj.frameInds = featuresCalc.setFrameInds(fullFrameNames, obj.frameNames);

            switch obj.feature
                case featuresEnums.cartVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dx(:, obj.frameInds));
                    
                case featuresEnums.cartAccelSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddx;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.ddx(:, obj.frameInds));
                    
                case featuresEnums.cartJerkSumSqu
                    obj.bf(1) = basisFeaturesEnums.dddx;
                    obj.bf(2) = basisFeaturesEnums.ddx;
                    obj.bf(3) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dddx(:, obj.frameInds));
                    
                case featuresEnums.cartCurvatureSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.bf(2) = basisFeaturesEnums.ddx;
                    obj.bf(3) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartCurvature(basisFeatures));
                    
                case featuresEnums.cartRadCurvatureSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.bf(2) = basisFeaturesEnums.ddx;
                    obj.bf(3) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartCurvature(basisFeatures).^(-1));

                case featuresEnums.cartBoundingBoxSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartBoundingBox(basisFeatures));
                    
                case featuresEnums.cartBoundingVolumeSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartBoundingVolume(basisFeatures));
                    
                case featuresEnums.cartDisplacementSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartDisplacement(basisFeatures));
                    if numel(obj.refFrameNames) ~= numel(obj.frameNames) 
                        % as in, there's many framenames but only one ref
                        for i = 2:length(obj.frameNames)
                            obj.refFrameNames{i} = obj.refFrameNames{1};
                        end
                    end
                    obj.referenceInds = featuresCalc.setFrameInds(fullFrameNames, obj.refFrameNames);
                    
                case featuresEnums.cartQuantityMotionSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartQuantityMotionSumSqu(basisFeatures));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.frameNames))*(1/numel(obj.frameNames));
                    end
                    
                case featuresEnums.cartWeightEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartWeightEffortSumSqu(basisFeatures));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.frameNames))*(1/numel(obj.frameNames));
                    end
                    if isempty(obj.windowLength) % equal weights
                        obj.windowLength = 1;
                    end
                    
                case featuresEnums.cartTimeEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartTimeFlowEffortSumSqu(basisFeatures.ddx));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.frameNames))*(1/numel(obj.frameNames));
                    end
                    if isempty(obj.windowLength) % equal weights
                        obj.windowLength = 1;
                    end
                    
                case featuresEnums.cartSpaceEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartSpaceEffortSumSqu(basisFeatures));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.frameNames))*(1/numel(obj.frameNames));
                    end
                    if isempty(obj.windowLength) % equal weights
                        obj.windowLength = 1;
                    end
                    
                case featuresEnums.cartFlowEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.dddx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartTimeFlowEffortSumSqu(basisFeatures.dddx));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.frameNames))*(1/numel(obj.frameNames));
                    end
                    if isempty(obj.windowLength) % equal weights
                        obj.windowLength = 1;
                    end
                    
                case featuresEnums.centreMassSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCentreOfMass(basisFeatures));
                    if isempty(obj.weights)
                        obj.weights = featuresCalc.calcMassWeights(obj.bodyNames, dynamicModel.model);
                    end
                    
                case featuresEnums.centreMassDisplacementSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCentreMassDisplacement(basisFeatures));
                    if isempty(obj.weights)
                        obj.weights = featuresCalc.calcMassWeights(obj.bodyNames, dynamicModel.model);
                    end
                    obj.referenceVals = featuresCalc.calcCentreMassReferencePose(fullFrameNames, obj.frameNames, dynamicModel); % determined by initial pose
                    
                case featuresEnums.angVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dq(:, obj.jointInds));
                    
                case featuresEnums.angAccelSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddq;
                    obj.bf(2) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.ddq(:, obj.jointInds));
                    
                case featuresEnums.angJerkSumSqu
                    obj.bf(1) = basisFeaturesEnums.dddq;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.bf(3) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngJerk(basisFeatures));
%                     obj.cf = @(basisFeatures) ...
%                         sumSqu(basisFeatures.dddq(:, obj.jointInds));       
                    
                case featuresEnums.angCurvatureSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngCurvature(basisFeatures));
                    
                case featuresEnums.angRadCurvatureSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngCurvature(basisFeatures).^(-1));
                    
                case featuresEnums.angQuantityMotionSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngQuantityMotionSumSqu(basisFeatures));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.jointNames))*(1/numel(obj.jointNames));
                    end
                   
                case featuresEnums.angWeightEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngWeightEffortSumSqu(basisFeatures));
                    if isempty(obj.weights) % equal weights
                        obj.weights = ones(size(obj.jointNames))*(1/numel(obj.jointNames));
                    end
                    if isempty(obj.windowLength) % equal weights
                        obj.windowLength = 1;
                    end
                    
                case featuresEnums.torqueSumSqu
                    obj.bf(1) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.tau(:, obj.jointInds));
                    
                case featuresEnums.torqueVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.dtau;
                    obj.bf(2) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dtau(:, obj.jointInds));
                    
                case featuresEnums.torqueAccelSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddtau;
                    obj.bf(2) = basisFeaturesEnums.dtau;
                    obj.bf(3) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.ddtau(:, obj.jointInds));
                    
                case featuresEnums.kineticEnergySumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.bf(2) = basisFeaturesEnums.M;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfKineticEnergy(basisFeatures));
                    
                case featuresEnums.angPowerSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.bf(2) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dq(:, obj.jointInds) .* basisFeatures.tau(:, obj.jointInds));
                    
                case featuresEnums.extensivenessMaxSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        maxAbs(obj.calcCfExtenstiveness(basisFeatures));
                    if isempty(obj.weights)
                        obj.weights = featuresCalc.calcMassWeights(obj.bodyNames, dynamicModel.model);
                    end
                    
                case featuresEnums.extensivenessSumSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfExtenstiveness(basisFeatures));
                    if isempty(obj.weights)
                        obj.weights = featuresCalc.calcMassWeights(obj.bodyNames, dynamicModel.model);
                    end
                    
                case featuresEnums.shapeDirectionSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        maxAbs(obj.calcCfSpaceDirectional(basisFeatures));
                    
                case featuresEnums.cartDistToTarget
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures)...
                        sumSqu(obj.calcCfCartDistToTarget(basisFeatures));
                    if numel(obj.refFrameNames) ~= numel(obj.frameNames) 
                        % as in, there's many framenames but only one ref
                        for i = 2:length(obj.frameNames)
                            obj.refFrameNames{i} = obj.refFrameNames{1}; % TO REVIEW
                        end
                    end
                    obj.referenceInds = featuresCalc.setFrameInds(fullFrameNames, obj.refFrameNames);
                    
                case featuresEnums.rotDistToTarget
                    obj.bf(1) = basisFeaturesEnums.R;
                    obj.cf = @(basisFeatures)...
                        sumSqu(obj.calcCfRotDistToTarget(basisFeatures));
                    if numel(obj.refFrameNames) ~= numel(obj.frameNames) 
                        % as in, there's many framenames but only one ref
                        for i = 2:length(obj.frameNames)
                            obj.refFrameNames{i} = obj.refFrameNames{1}; % TO REVIEW
                        end
                    end
                    obj.referenceInds = featuresCalc.setFrameInds(fullFrameNames, obj.refFrameNames);
                              
                otherwise
                    error('featuresCalc not defined');
            end
        end
        
        function specialize(obj, fullJointNames, fullFrameNames, dynamicModel, trajU, trajX)
            switch obj.feature
                case featuresEnums.cartDistToTarget
                    % Get cartesian position of reference frame at given timeframe
                    lenDofs = lenght(obj.refFrameNames);
                    lenTimes = lenght(obj.referenceTimes);
                    obj.referenceVals = zeros(lenTimes, lenDofs * 3);
                    default_state = dynamicModel.getState();
                    
                    for i = 1:lenTimes
                        time = obj.referenceTimes(i);
                        targetState = trajX(time,:);
                        dynamicModel.updateState(targetState(1:size(targetState,1)/2),targetState(size(targetState,1)/2+1:end));
                        dynamicModel.forwardKinematics();
                        
                        for j=1:lenDofs
                            inds = obj.referenceInds((1:3)+(j-1)*3);
                            targetInd = find(ismember(fullFrameNames, obj.refFrameNames(j)));
                            obj.referenceVals(i,inds) = dynamicModel.getEndEffectorPosition(targetInd);
                        end
                    end
                    % Set dynamic model state to default
                    dynamicModel.updateState(default_state(1:length(default_state)/2),default_state(length(default_state)/2+1:end));
                    
               case featuresEnums.rotDistToTarget
                    % Get cartesian position of reference frame at given timeframe
                    lenDofs = lenght(obj.refFrameNames);
                    lenTimes = lenght(obj.referenceTimes);
                    obj.referenceVals = zeros(3, lenDofs * 3, lenTimes);
                    default_state = dynamicModel.getState();
                    
                    for i = 1:lenTimes
                        time = obj.referenceTimes(i);
                        targetState = trajX(time,:);
                        dynamicModel.updateState(targetState(1:size(targetState,1)/2),targetState(size(targetState,1)/2+1:end));
                        dynamicModel.forwardKinematics();
                            for j=1:lenDofs
                                inds = obj.referenceInds((1:3)+(j-1)*3);
                                targetInd = find(ismember(fullFrameNames, obj.refFrameNames(j)));
                                obj.referenceVals(:,inds, i) = dynamicModel.getEndEffectorPosition(targetInd);
                            end
                    end
                    % Set dynamic model state to default
                    dynamicModel.updateState(default_state(1:length(default_state)/2),default_state(length(default_state)/2+1:end)); 
                
                otherwise
                    error('featuresCalc not defined');
            end
        end
        
        
        function cf = calcFeature(obj, features)
            % use the features calculated from the
            cf = obj.cf(features) / obj.normCoeff;
        end
        
        function cf = calcCfCartCurvature(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, lenDofs);
            for i = 1:lenTime
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    num = norm(cross(features.ddx(i, inds), features.dx(i, inds))); % todo can vectorize
                    den = norm(features.dx(i, inds).^(3/2));
                    if den > 0
                        cfVal = num/den; % if there is no dx, then we're dividing by zero
                    else
                        cfVal = 0;
                    end
                    
                    cf(i, j) = cfVal;
                end
            end
        end
        
        function cf = calcCfAngJerk(obj, features)
            cf = features.dddq(:, obj.jointInds);       
        end
        
        function cf = calcCfAngCurvature(obj, features)
            cf = features.ddq ./ (features.dq.^2); % using the general form of the curvature formula
            divByZero = isinf(cf);
            cf(divByZero) = 0;
        end

        function cf = calcCfCartBoundingBox(obj, features)
            lenTime = size(features.x, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 3);
            for i = 1:lenTime
                allCartPos = zeros(lenDofs, 3);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    allCartPos(j, :) = features.x(i, inds);
                end
                
                cfVal = zeros(1, 3);
                for j = 1:3
                    cartRange = max(allCartPos(:, j)) - min(allCartPos(:, j));
                    cfVal(j) = cartRange;
                end
                
                cf(i, :) = cfVal;
            end
        end

        
        function cf = calcCfCartBoundingVolume(obj, features)
            cfValTemp = calcCfCartBoundingBox(obj, features);
            cf = cfValTemp(:, 1).*cfValTemp(:, 2).*cfValTemp(:, 3);
        end
        
        function cf = calcCfCartDisplacement(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, lenDofs);
            for i = 1:lenTime
                for j = 1:lenDofs
                    checkInds = obj.frameInds((1:3)+(j-1)*3);
                    refInds = obj.referenceInds((1:3)+(j-1)*3);
                    cfVal = norm(features.x(i, checkInds) - features.x(i, refInds));
                    cf(i, j) = cfVal;
                end
            end
        end
        
        function cf = calcCfCentreOfMass(obj, features)
            lenTime = size(features.x, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 3);
            for i = 1:lenTime
                cfVal = zeros(1, 3);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    cfVal = cfVal + obj.weights(j)*features.x(i, inds);
                end
                
                com = cfVal / sum(obj.weights);
                cf(i, :) = com;
            end
        end
        
        function cf = calcCfCartQuantityMotionSumSqu(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 3);
            for i = 1:lenTime
                cfVal = zeros(1, 3);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    cfVal = cfVal + obj.weights(j)*features.dx(i, inds);
                end
                
                cfNorm = cfVal / sum(obj.weights);
                cf(i, :) = cfNorm;
            end
        end
        
        function cf = calcCfAngQuantityMotionSumSqu(obj, features)
            lenTime = size(features.dq, 1);
            lenDofs = size(obj.jointInds, 2);
            
            cf = zeros(lenTime, lenDofs);
            for i = 1:lenTime
                cfVal = 0;
                for j = 1:lenDofs
                    cfVal = cfVal + obj.weights(j)*features.dq(i, j);
                end
                
                com = cfVal / sum(obj.weights);
                cf(i, :) = com;
            end
        end
        
        function cf = calcCfCartWeightEffortSumSqu(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cfCum = zeros(lenTime, 3);
            cf = zeros(lenTime, 1);
           
            for i = 1:lenTime  % calc the sum over all the joints
                cfVal = zeros(1, 3);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    cfVal = cfVal + obj.weights(j)*norm(features.dx(i, inds).^ 2);
                end
                
                cfCum(i, :) = cfVal;
            end
            
            for i = 1:lenTime   % then find the max of a sliding window
                indMax = (i-obj.windowLength+1):i; % that is, if winlen = 1, then indMax is just i
                if min(indMax) < 1
                    indMax = 1:i;
                end
                
                valsToCheck = vecnorm(cfCum(indMax, :), 2, 2);
                cf(i) = max(valsToCheck); 
            end
        end    
        
        function cf = calcCfAngWeightEffortSumSqu(obj, features)
            lenTime = size(features.dq, 1);
            lenDofs = size(obj.jointInds, 2);
            
            cfCum = zeros(lenTime, 1);
            cf = zeros(lenTime, 1);
            for i = 1:lenTime
                cfVal = 0;
                for j = 1:lenDofs
                    cfVal = cfVal + obj.weights(j)*norm(features.dq(i, j).^ 2);
                end
                
                cfCum(i, :) = cfVal;
            end
            
            for i = 1:lenTime
                indMax = (i-obj.windowLength+1):i; % that is, if winlen = 1, then indMax is just i
                if min(indMax) < 1
                    indMax = 1:i;
                end
                
                valsToCheck = vecnorm(cfCum(indMax, :), 2, 2);
                cf(i) = max(valsToCheck); 
            end
        end
        
%         function cf = calcCfCartTimeEffortSumSqu(obj, features)
%             lenTime = size(features.ddx, 1);
%             lenDofs = size(obj.frameInds, 2)/3;
%             
%             cf = zeros(lenTime, 1);
%             for i = 1:lenTime
%                 cfVal = zeros(1, 3);
%                 for j = 1:lenDofs % sum over all the joints at time i
%                     inds = obj.frameInds((1:3)+(j-1)*3);
%                     cfVal = cfVal + obj.weights(j)*features.ddx(i, inds);
%                 end
%                 
%                 cfCum(i, :) = cfVal; 
%             end
%             
%             for i = 1:lenTime
%                 indMax = (i-obj.windowLength+1):i; % that is, if winlen = 1, then indMax is just i
%                 if min(indMax) < 1
%                     indMax = 1:i;
%                 end
%                 
%                 valsToCheck = vecnorm(cfCum(indMax, :), 2, 2);
%                 cf(i) = sum(valsToCheck)/obj.windowLength;
%             end
%         end
        
        function cf = calcCfCartSpaceEffortSumSqu(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cfVal = zeros(lenTime, lenDofs);
            
            for i = (obj.windowLength):lenTime
                for k = 1:lenDofs % the dofs
                    inds = obj.frameInds((1:3)+(k-1)*3);
                    
                    if obj.windowLength == 1
                        currCf = 1;
                    else
                        j_start = i - obj.windowLength+1;
                        j_end = i;
                        space_k = 0;
                        
                        if j_start < 0
                            j_start = 0;
                        end
                        
                        for j = (j_start+1):j_end % 1:T
                            space_k = space_k + norm(features.x(j, inds) - features.x(j-1, inds));
                        end
                        
                        delta = norm(features.x(j_end, inds) - features.x(j_start, inds));
                        currCf = obj.weights(k)*space_k/delta;
                    end
                    
                    
                    cfVal(i, k) = currCf;
                end
            end
            
            cf = sum(cfVal, 2);
        end
        
        function cf = calcCfCartTimeFlowEffortSumSqu(obj, features_x) % time and flow effort are calculated the same way
            lenTime = size(features_x, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 1);
            for i = 1:lenTime
                cfVal = zeros(1, 3);
                for j = 1:lenDofs % sum over all the joints at time i
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    cfVal = cfVal + obj.weights(j)*features_x(i, inds);
                end
                
                cfCum(i, :) = cfVal;
            end
            
            for i = 1:lenTime
                indMax = (i-obj.windowLength+1):i; % that is, if winlen = 1, then indMax is just i
                if min(indMax) < 1
                    indMax = 1:i;
                end
                
                valsToCheck = vecnorm(cfCum(indMax, :), 2, 2);
                cf(i) = sum(valsToCheck)/obj.windowLength;
            end
        end
        
        function cf = calcCfCentreMassDisplacement(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 3);
            for i = 1:lenTime
                cfVal = zeros(3, 1);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    cfVal = features.x(i, inds) - obj.referenceVals(inds);
                end
                
                cf(i, :) = cfVal;
            end
        end
        
        function cf = calcCfShapeDirectionSumSqu(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, lenDofs);
            for i = 2:lenTime
                cfVal = zeros(3, 1);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    temp(i) = norm(obj.weights(j)*features.x(i, inds) - obj.weights(j)*features.x(i-1, inds));
                end
                
                cf(lenTime, j) = obj.weight(j)*sum(temp);
                
                cf(i, j) = cfVal;
            end
        end
        
        function cf = calcCfExtenstiveness(obj, features)
            lenTime = size(features.x, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            com = obj.calcCfCentreOfMass(features);
            
            cf = zeros(lenTime, lenDofs);
            for i = 1:lenTime
                cfVal = zeros(3, 1);
                for j = 1:lenDofs
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    cf(i, j) = obj.weights(j)*norm(com(i, :) - features.x(i, inds));
                end
            end
        end
        
        function cf = calcCfSpaceDirectional(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, lenDofs);
            for i = 2:lenTime
                for j = 1:lenDofs
                    % assemble
                    inds = obj.frameInds((1:3)+(j-1)*3);
                    tempx(j, :) = features.x(i, inds);
                    tempdx(j, :) = features.dx(i, inds);
                    tempddx(j, :) = features.ddx(i, inds);
                end
                
                coeff = pca(tempx)';
                
                pcadx = zeros(lenDofs, 3);
                pcaddx = zeros(lenDofs, 3);
                for j = 1:lenDofs
                    pcadx(j, :) = coeff.*tempdx(j, :);
                    pcaddx(j, :) = coeff.*tempddx(j, :);
                end
           
                num = (pcaddx(2)*pcadx(1) - pcaddx(1)*pcadx(2)) .^ 2;
                den = ((pcadx(1) + pcadx(2)) .^2) .^ (3/2);
                
                if den ~= 0
                    cf(i, j) = num/den;
                else
                    cf(i, j) = 0;
                end
            end
        end
        
        function cf = calcCfKineticEnergy(obj, features)
            lenTime = size(features.dq, 1);
            cf = zeros(lenTime, 1);
            
            for i = 1:lenTime
                cf(i, :) = (features.dq(i, :) * features.M(:, :, i) * features.dq(i, :)')';
            end
        end
        
        function cf = calcCfCartDisToTarget(obj, features)
            lenTime = size(features.x, 1);
            lenRefTimes = length(obj.referenceTimes); % We assume refFrameNames and frameNames have same lenght and share same order
            lenDofs = length(obj.refFrameNames);
            cf = zeros(lenTime, 1);
            
            for i = 1:lenTime
                cfVal = 0;
                for j = 1:lenRefTimes
                    for k= 1:lenDofs
                        inds = obj.frameInds((1:3)+(k-1)*3);
                        rInds = obj.referenceInds((1:3)+(k-1)*3);
                        curPose = features.x(i,inds);
                        refPose = obj.refVals(j,rInds);
                        % Here we consider that distances are computed by
                        % tuples (i.e., 1st frame in frameInds with 1st frame
                        % in refFrames). In case of several target poses, we
                        % add all distances (this definition needs to be
                        % reviewed)
                        cfVal = cfVal + norm(refPose - curPose); 
                    end
                end
                cf(i,:) = cfVal;
            end
        end
        
        function cf = calcCfRotDisToTarget(obj, features)
            lenTime = size(features.x, 1);
            lenRefTimes = length(obj.referenceTimes); % We assume refFrameNames and frameNames have same lenght and share same order
            lenDofs = length(obj.refFrameNames);
            cf = zeros(lenTime, 1);
            
            for i = 1:lenTime
                cfVal = 0;
                for j = 1:lenRefTimes
                    for k= 1:lenDofs
                        inds = obj.frameInds((1:3)+(k-1)*3);
                        rInds = obj.referenceInds((1:3)+(k-1)*3);
                        curPose = features.R(:,inds,i);
                        refPose = obj.refVals(:,rInds,j);
                        cfVal = cfVal + acos((trace(curPose*refPose') -1)/2); 
                    end
                end
                cf(i,:) = cfVal;
            end
        end
    end
    
    methods(Static=true)
        function jointInds = setJointInds(fullJointNames, jointNames)
            jointInds = [];
            for i = 1:length(jointNames)
                jointInds(i) = find(ismember(fullJointNames, jointNames(i)));
            end
        end
        
        function frameInds = setFrameInds(fullFrameNames, frameNames)
            frameInds = [];
            for i = 1:length(frameNames)
                inds = (1:3)+(i-1)*3;
                targetInd = find(ismember(fullFrameNames, frameNames(i)));
                frameInds(inds) = (1:3)+(targetInd-1)*3;
            end
        end
        
        function weights = calcMassWeights(bodyNames, rlModel)
            weights = [];
            allBodyNames = {rlModel.bodies.name};
            for i = 1:length(bodyNames)
                findInd = find(ismember(allBodyNames, bodyNames{i}));
                weights = [weights rlModel.bodies(findInd).m];
            end
            
            weights = weights / sum(weights);
        end
        
        function x = calcCentreMassReferencePose(frameNamesFull, frameNames, dynamicModel)
            for j = 1:length(frameNames)
                inds = (1:3)+(j-1)*3;
                
                sourceInds = find(ismember(frameNamesFull, frameNames{j}));
                x(inds) = dynamicModel.getEndEffectorPosition(sourceInds);
            end
        end
                
%         function frameInds = parseReferenceXYZ(referenceCell)
%             frameInds = [];
%             if ismember(referenceCell, "x")
%                 frameInds = [frameInds 1];
%             end
%             
%             if ismember(referenceCell, "y")
%                 frameInds = [frameInds 2];
%             end
%             
%             if ismember(referenceCell, "z")
%                 frameInds = [frameInds 3];
%             end
%         end
    end
end