classdef featuresCalc < handle
    % denotes the type of features to calculate from the dynamic model.
    % this class is used in conjunction with IOCInstance
    properties
        % friendly name for the feature
        name
        % loaded feature type from .json. this is a enum of featureEnums
        feature
        % loaded parameter form .json. these are modifiers that are parsed 
        % by this function and IOCInstance to determine how to calculate 
        % the cost function
        
        jointNames % which joint frame to use. leave blank if want to use all the
                   % joint frames
        frameNames % which transform frame to use. 
                   % for COM-related functions, leave blank if want to use all the body frames
%         bodyNames  
        refFrameNames % reference frame or pose. its use depends on the cf
        weights % weights for cfs that require it
        axes % xyz for bounding box or other directional-able features
        
        % determines if the frameInds belongs to frames (ie cartesian) or
        % joints (ie revolute)
        jointType
        % frameInds of the full array used. this is calculated by
        % IOCInstance from obj.param
        jointInds
        frameInds
%         bodyInds
        refInds
        refVals
        refTimes
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
            addOptional(externalParamParser, 'name', []);
            addOptional(externalParamParser, 'feature', []);
            addOptional(externalParamParser, 'jointNames', []);
            addOptional(externalParamParser, 'frameNames', []);
%             addOptional(externalParamParser, 'bodyNames', []);
            addOptional(externalParamParser, 'weights', []);
            addOptional(externalParamParser, 'axes', []);
            addOptional(externalParamParser, 'refFrameNames', []);            
            addOptional(externalParamParser, 'windowLength', []);
            addOptional(externalParamParser, 'normCoeff', 1);
            addOptional(externalParamParser, 'refTimes', []);

            parse(externalParamParser, jsonBlob);
            jsonParsed = externalParamParser.Results;
            
            obj.name = jsonParsed.name;
            if isempty(obj.name)
                obj.name = jsonParsed.feature; % give the function a default name
            end
            
            obj.feature = featuresEnums.(jsonParsed.feature); % enum parse
            obj.jointNames = jsonParsed.jointNames(:);
            obj.frameNames = jsonParsed.frameNames(:);
%             obj.bodyNames = jsonParsed.bodyNames(:);
            obj.refFrameNames = jsonParsed.refFrameNames;
            obj.refTimes = jsonParsed.refTimes(:);
            
            obj.weights = jsonParsed.weights(:)';
            obj.axes = jsonParsed.axes(:)';
            obj.windowLength = jsonParsed.windowLength;
            
            obj.normCoeff = jsonParsed.normCoeff;
            
            if numel(obj.refFrameNames) ~= numel(obj.frameNames)
                % as in, there's many framenames but only one ref
                for i = 2:length(obj.frameNames)
                    obj.refFrameNames{i} = obj.refFrameNames{1};
                end
            end
        end
        
        function initCf(obj, dynamicModel, fullJointNames)
            sumSqu = @(x) sum(x.^2, 2);
%             sumAbs = @(x) sum(abs(x), 2);
%             maxAbs = @(x) max(abs(x), [], 2);
            
            obj.bf = basisFeaturesEnums.dddx;  % initializing to something temp for now

            switch obj.feature
                case featuresEnums.cartVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dx(:, obj.frameInds) .* repmat(obj.axes, size(basisFeatures.dx, 1), 1));
                    obj.checkDefaultAxes();
                    
                case featuresEnums.cartAccelSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddx;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.ddx(:, obj.frameInds) .* repmat(obj.axes, size(basisFeatures.ddx, 1), 1));
                    obj.checkDefaultAxes();
                    
                case featuresEnums.cartJerkSumSqu
                    obj.bf(1) = basisFeaturesEnums.dddx;
                    obj.bf(2) = basisFeaturesEnums.ddx;
                    obj.bf(3) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dddx(:, obj.frameInds) .* repmat(obj.axes, size(basisFeatures.dddx, 1), 1));
                    obj.checkDefaultAxes();
                    
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
                    
                case featuresEnums.cartQuantityMotionSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfWeightedEndEffVelo(basisFeatures));
                    obj.checkDefaultWeights(obj.frameNames);
                    
                case featuresEnums.cartWeightEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartWeightEffortSumSqu(basisFeatures));
                    obj.checkDefaultWeights(obj.frameNames);
                    obj.checkDefaultLengths();
                    
                case featuresEnums.cartTimeEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartTimeFlowEffortSumSqu(basisFeatures.ddx));
                    obj.checkDefaultWeights(obj.frameNames);
                    obj.checkDefaultLengths();
                    
                case featuresEnums.cartSpaceEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartSpaceEffortSumSqu(basisFeatures));
                    obj.checkDefaultWeights(obj.frameNames);
                    obj.checkDefaultLengths();
                    
                case featuresEnums.cartFlowEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.dddx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCartTimeFlowEffortSumSqu(basisFeatures.dddx));
                    obj.checkDefaultWeights(obj.frameNames);
                    obj.checkDefaultLengths();
                    
                case featuresEnums.centreMassSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfWeightedEndEffPos(basisFeatures));
                    obj.checkDefaultWeightsCOM(dynamicModel);
                    obj.checkDefaultAxes();
                    
                case featuresEnums.centreMassVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfWeightedEndEffVelo(basisFeatures));
                    obj.checkDefaultWeightsCOM(dynamicModel);
                    obj.checkDefaultAxes();                    
                    
                case featuresEnums.centreMassDisplacementSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCentreMassDisplacement(basisFeatures));
                    obj.checkDefaultWeightsCOM(dynamicModel);
                    obj.refVals = featuresCalc.calcCentreMassReferencePose(fullFrameNames, obj.frameNames, dynamicModel); % determined by initial pose
                    
                case featuresEnums.angVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dq(:, obj.jointInds));
                    obj.checkDefaultJoints(fullJointNames);
                    
                case featuresEnums.angAccelSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddq;
                    obj.bf(2) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.ddq(:, obj.jointInds));
                    obj.checkDefaultJoints(fullJointNames);
                    
                case featuresEnums.angJerkSumSqu
                    obj.bf(1) = basisFeaturesEnums.dddq;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.bf(3) = basisFeaturesEnums.dq;
%                     obj.cf = @(basisFeatures) ...
%                         sumSqu(obj.calcCfAngJerk(basisFeatures));
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dddq(:, obj.jointInds));       
                    obj.checkDefaultJoints(fullJointNames);
                    
                case featuresEnums.angCurvatureSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngCurvature(basisFeatures));
                    obj.checkDefaultJoints(fullJointNames);
                    
                case featuresEnums.angRadCurvatureSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.bf(2) = basisFeaturesEnums.ddq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngCurvature(basisFeatures).^(-1));
                    obj.checkDefaultJoints(fullJointNames);
                    
                case featuresEnums.angQuantityMotionSumSqu
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngQuantityMotionSumSqu(basisFeatures));
                    obj.checkDefaultWeights(obj.frameNames);
                   
                case featuresEnums.angWeightEffortSumSqu
                    obj.bf(1) = basisFeaturesEnums.dq;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfAngWeightEffortSumSqu(basisFeatures));
                    obj.checkDefaultWeights(obj.jointNames);
                    obj.checkDefaultLengths();
                    
                case featuresEnums.torqueSumSqu
                    obj.bf(1) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.tau(:, obj.jointInds));
                    obj.checkDefaultJoints();
                    
                case featuresEnums.torqueVeloSumSqu
                    obj.bf(1) = basisFeaturesEnums.dtau;
                    obj.bf(2) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.dtau(:, obj.jointInds));
                    obj.checkDefaultJoints();
                    
                case featuresEnums.torqueAccelSumSqu
                    obj.bf(1) = basisFeaturesEnums.ddtau;
                    obj.bf(2) = basisFeaturesEnums.dtau;
                    obj.bf(3) = basisFeaturesEnums.tau;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(basisFeatures.ddtau(:, obj.jointInds));
                    obj.checkDefaultJoints();
                    
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
                    obj.checkDefaultJoints();
                    
                case featuresEnums.extensivenessMaxSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        maxAbs(obj.calcCfExtenstiveness(basisFeatures));
                    obj.checkDefaultWeightsCOM(dynamicModel);
                    
                case featuresEnums.extensivenessSumSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfExtenstiveness(basisFeatures));
                    obj.checkDefaultWeightsCOM(dynamicModel);
                    
                case featuresEnums.shapeDirectionSumSqu
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures) ...
                        maxAbs(obj.calcCfSpaceDirectional(basisFeatures));
                    
                case featuresEnums.cartDistToTarget
                    obj.bf(1) = basisFeaturesEnums.x;
                    obj.cf = @(basisFeatures)...
                        sumSqu(obj.calcCfCartDistToTarget(basisFeatures));
                    
                case featuresEnums.rotDistToTarget
                    obj.bf(1) = basisFeaturesEnums.R;
                    obj.cf = @(basisFeatures)...
                        sumSqu(obj.calcCfRotDistToTarget(basisFeatures));
                               
                case featuresEnums.centreMassVeloRelativeToFrame
                    obj.bf(1) = basisFeaturesEnums.dx;
                    obj.cf = @(basisFeatures) ...
                        sumSqu(obj.calcCfCentreMassDisplacementRelativeToFrame(basisFeatures));
                    obj.checkDefaultWeightsCOM(dynamicModel);
                    obj.checkDefaultAxes();
                    
                otherwise
                    error('featuresCalc not defined');
            end
        end
        
        function initInds(obj, fullJointNames, fullFrameNames)  
            obj.jointInds = featuresCalc.setJointInds(fullJointNames, obj.jointNames);
            obj.frameInds = featuresCalc.setFrameInds(fullFrameNames, obj.frameNames); 
            
            obj.refInds = featuresCalc.setFrameInds(fullFrameNames, obj.refFrameNames);
        end
        
        function checkDefaultWeights(obj, frameJointNames)
            if isempty(obj.weights) % equal weights
                fprintf('featuresCalc: %s (%s): No weights defined in json, using uniform weights over all frames/joints\n', obj.name, obj.feature);
                obj.weights = ones(size(frameJointNames))*(1/numel(frameJointNames));
            end
        end
        
        function checkDefaultWeightsCOM(obj, dynamicModel)
            if isempty(obj.frameNames)
                fprintf('featuresCalc: %s (%s): No frames defined in json, using mass of individual links\n', obj.name, obj.feature);
                obj.frameNames = {dynamicModel.model.bodies.name}';
            else
                fprintf('featuresCalc: %s (%s): Frames are defined in json, verify that frame was intentional and not for refFrame\n', obj.name, obj.feature);
            end
            
            if isempty(obj.weights)
                fprintf('featuresCalc: %s (%s): No weights defined in json, using mass of all listed bodynames\n', obj.name, obj.feature);
                obj.weights = featuresCalc.calcMassWeights(obj.frameNames, dynamicModel.model);
            end
            
            % if an entry does not have any weights, remove it from the
            % frame names to reduce calculation costs
            findFilledInds = find(obj.weights > 0);
            obj.weights = obj.weights(findFilledInds);
            obj.frameNames = obj.frameNames(findFilledInds);
        end
        
        function checkDefaultAxes(obj)
            if isempty(obj.axes)
                fprintf('featuresCalc: %s (%s): No axes defined in json, using [1 1 1]\n', obj.name, obj.feature);
                obj.axes = ones(1, 3);
            end
        end
        
        function checkDefaultLengths(obj)
            if isempty(obj.windowLength) % equal weights
                fprintf('featuresCalc: %s (%s): No length defined in json, using 1\n', obj.name, obj.feature);
                obj.windowLength = 1;
            end
        end
        
        function checkDefaultJoints(obj, fullJointNames)
            if isempty(obj.jointNames) % if no joint names defined, init to all the joints
                fprintf('featuresCalc: %s (%s): No joints defined in json, using all joints\n', obj.name, obj.feature);
                obj.jointNames = fullJointNames;
            end
        end
        
        function specialize(obj, fullJointNames, fullFrameNames, dynamicModel, trajU, trajX)
            switch obj.feature
                case featuresEnums.cartDistToTarget
                    % Get cartesian position of reference frame at given timeframe
                    lenDofs = lenght(obj.refFrameNames);
                    lenTimes = lenght(obj.refTimes);
                    obj.refVals = zeros(lenTimes, lenDofs * 3);
                    default_state = dynamicModel.getState();
                    
                    for i = 1:lenTimes
                        time = obj.refTimes(i);
                        targetState = trajX(time,:);
                        dynamicModel.updateState(targetState(1:size(targetState,1)/2),targetState(size(targetState,1)/2+1:end));
                        dynamicModel.forwardKinematics();
                        
                        for j=1:lenDofs
                            inds = obj.refInds(j, :);
                            targetInd = find(ismember(fullFrameNames, obj.refFrameNames(j)));
                            obj.refVals(i,inds) = dynamicModel.getEndEffectorPosition(targetInd);
                        end
                    end
                    % Set dynamic model state to default
                    dynamicModel.updateState(default_state(1:length(default_state)/2),default_state(length(default_state)/2+1:end));
                    
               case featuresEnums.rotDistToTarget
                    % Get cartesian position of reference frame at given timeframe
                    lenDofs = lenght(obj.refFrameNames);
                    lenTimes = lenght(obj.refTimes);
                    obj.refVals = zeros(3, lenDofs * 3, lenTimes);
                    default_state = dynamicModel.getState();
                    
                    for i = 1:lenTimes
                        time = obj.refTimes(i);
                        targetState = trajX(time,:);
                        dynamicModel.updateState(targetState(1:size(targetState,1)/2),targetState(size(targetState,1)/2+1:end));
                        dynamicModel.forwardKinematics();
                            for j=1:lenDofs
                                inds = obj.refInds(j, :);
                                targetInd = find(ismember(fullFrameNames, obj.refFrameNames(j)));
                                obj.refVals(:,inds, i) = dynamicModel.getEndEffectorPosition(targetInd);
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
                    inds = obj.frameInds(j, :);
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
                    inds = obj.frameInds(j, :);
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
                    checkInds = obj.frameInds(j, :);
                    refInds = obj.refInds(j, :);
                    cfVal = norm(features.x(i, checkInds) - features.x(i, refInds));
                    cf(i, j) = cfVal;
                end
            end
        end
        
        function cf = calcCfCentreMassDisplacementRelativeToFrame(obj, features)
            lenTime = size(features.dx, 1);
            cf = zeros(lenTime, 3);

            cf_com = obj.calcCfWeightedEndEffPos(features);
            for i = 1:lenTime
                cfNew = cf_com(i, :) - features.x(i, obj.refInds);
                cf(i, :) = cfNew.*obj.axes;
            end
        end
        
        function cf = calcCfWeightedEndEffPos(obj, features)
            lenTime = size(features.x, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 3);
            for i = 1:lenTime
                cfVal = zeros(1, 3);
                for j = 1:lenDofs
                    inds = obj.frameInds(j, :);
                    cfVal = cfVal + obj.weights(j)*features.x(i, inds);
                end
                
                com = cfVal / sum(obj.weights);
                cf(i, :) = com.*obj.axes;
            end
        end
        
        function cf = calcCfWeightedEndEffVelo(obj, features)
            lenTime = size(features.dx, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            cf = zeros(lenTime, 3);
            for i = 1:lenTime
                cfVal = zeros(1, 3);
                for j = 1:lenDofs
                    inds = obj.frameInds(j, :);
                    cfVal = cfVal + obj.weights(j)*features.dx(i, inds);
                end
                
                cfNorm = cfVal / sum(obj.weights);
                cf(i, :) = cfNorm .* obj.axes;
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
                    inds = obj.frameInds(j, :);
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
%                     inds = obj.frameInds(j, :);
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
                    inds = obj.frameInds(j, :);
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
                    inds = obj.frameInds(j, :);
                    cfVal = features.x(i, inds) - obj.refVals(inds);
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
                    inds = obj.frameInds(j, :);
                    temp(i) = norm(obj.weights(j)*features.x(i, inds) - obj.weights(j)*features.x(i-1, inds));
                end
                
                cf(lenTime, j) = obj.weight(j)*sum(temp);
                
                cf(i, j) = cfVal;
            end
        end
        
        function cf = calcCfExtenstiveness(obj, features)
            lenTime = size(features.x, 1);
            lenDofs = size(obj.frameInds, 2)/3;
            
            com = obj.calcCfWeightedEndEffPos(features);
            
            cf = zeros(lenTime, lenDofs);
            for i = 1:lenTime
                cfVal = zeros(3, 1);
                for j = 1:lenDofs
                    inds = obj.frameInds(j, :);
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
                    inds = obj.frameInds(j, :);
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
            lenRefTimes = length(obj.refTimes); % We assume refFrameNames and frameNames have same lenght and share same order
            lenDofs = length(obj.refFrameNames);
            cf = zeros(lenTime, 1);
            
            for i = 1:lenTime
                cfVal = 0;
                for j = 1:lenRefTimes
                    for k= 1:lenDofs
                        inds = obj.frameInds((1:3)+(k-1)*3);
                        rInds = obj.refInds((1:3)+(k-1)*3);
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
            lenRefTimes = length(obj.refTimes); % We assume refFrameNames and frameNames have same lenght and share same order
            lenDofs = length(obj.refFrameNames);
            cf = zeros(lenTime, 1);
            
            for i = 1:lenTime
                cfVal = 0;
                for j = 1:lenRefTimes
                    for k= 1:lenDofs
                        inds = obj.frameInds((1:3)+(k-1)*3);
                        rInds = obj.refInds((1:3)+(k-1)*3);
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
                targetInd = find(ismember(fullFrameNames, frameNames(i)));
%                 inds = (1:3)+(i-1)*3;
                frameInds(i,:) = (1:3)+(targetInd-1)*3;
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