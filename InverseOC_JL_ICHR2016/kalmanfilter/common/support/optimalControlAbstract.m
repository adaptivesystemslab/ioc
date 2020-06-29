classdef optimalControlAbstract < handle
    properties
        setupMode = '';
        
        c_cost = []; % fixed in DOC, changes in IOC
        c_const = [];
        featureSet = []; % fixed in IOC, changes in DOC
        normalizeCostSet = [];
        
        costFunctionStruct = [];

        splineDt;
        windowTime;
        splineWindowAdvance = 0;
        
        constPoints = [];
%         splinePoints = [];
        
        matVecStructDoc = [];
        matVecStructIoc;
        
        rmseFct = @(a,b) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2)));
    end
    
%     methods(Abstract = true)
%         J = calcCostFunctionIndividual(costFunctionName); % normalize features and calculate cost
%     end
    
    methods
        function obj = init(obj, setupMode, featureSet)
            obj.setupMode = setupMode;
            obj.featureSet = featureSet;
              
            % setup cost functions and normalization
            obj.setupCostFunctions();
            obj.setupCostNormalization();
            
            % setup time information
            obj.setupTimeArray();
            
            % setup constraint information
            obj.setupConstraintPoints();
            
            % if the feature set is empty, initialize the feature set based on
            % the constraint points
            if isempty(obj.featureSet.q)
                obj.featureSet.splineModelInit();
                obj.featureSet.splineModelSetViaPointsStruct(obj.constPoints);
                obj.featureSet.calcQFromSpline(obj.windowTime);
                
                obj.featureSet.dt = obj.splineDt;
            end
            
            % setup splining
            obj.featureSet.splineModelInit();
            splineViaPoints = obj.setupSplineViaPoints();
            obj.featureSet.splineModelSetViaPointsStruct(splineViaPoints);
        end
        
        function obj = docInit(obj)
%             obj.featureSet.time = obj.windowTime;
%             obj.featureSet.dt = obj.splineDt;
        end
        
        function obj = iocInit(obj)

        end
        
%         function obj = loadDataset(obj, featureSet, normalizeCostSet)
%             obj.featureSet = featureSet;
%             obj.normalizeCostSet = normalizeCostSet;
%         end
        
%         function obj = 
        
        function [J_cost_sum, J_cost_array] = calcCostFunction(obj, c_cost, c_const)
            J_cost_array = zeros(size(obj.costFunctionStruct));
            for i = 1:length(obj.costFunctionStruct)
                J_cost_array(i) = obj.costFunctionStruct(i).function(obj.featureSet, obj.normalizeCostSet);
            end
            
            J_cost_array_weighted = c_cost .* J_cost_array;
            J_cost_sum = sum(J_cost_array_weighted);
        end
        
        function [ceq_vec, ceq_q, ceq_dq, ceq_ddq] = calcEqualityConstFunction(obj, c_cost, c_const)
            % assuming that the constraints consist of q, dq, ddq
            constInds = [obj.constPoints.inds];
            
            constArray = [obj.constPoints.q];
            featureArray = obj.featureSet.q(constInds, :);
            ceq_q = c_const(1)*(featureArray - constArray);
            
            constArray = [obj.constPoints.dq];
            featureArray = obj.featureSet.dq(constInds, :);
            ceq_dq = c_const(2)*(featureArray - constArray);

            constArray = [obj.constPoints.ddq];
            featureArray = obj.featureSet.ddq(constInds, :);
            ceq_ddq = c_const(3)*(featureArray - constArray);
            
            ceq_mat = [ceq_q ceq_dq ceq_ddq];
            ceq_vec = mat2vec(ceq_mat);
        end
        
        function obj = yVecUpdateFeatureSet(obj, y_opt_vec)
            y_opt_mat = vec2mat_multi(y_opt_vec, obj.matVecStructDoc); % convert the vector back into matrix form for reading
            obj.featureSet.splineModelUpdateViaPoints(y_opt_mat);
            obj.featureSet.calcQFromSpline(obj.windowTime);
            obj.featureSet.calcFeaturesFromQ();
        end
        
        function optStruct = getIocOptStruct(obj)
            optStruct = obj.iocOptNormStruct;
        end
        
        function optStruct = getIocOptStructBestPivot(obj)
            optStruct = obj.iocOptNormStruct(obj.optPivotInd);
        end
        
        function optStruct = getDocOptStruct(obj)
            optStruct = obj.docOptStruct;
        end
    end
end