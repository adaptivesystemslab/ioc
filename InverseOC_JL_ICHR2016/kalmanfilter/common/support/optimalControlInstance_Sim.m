classdef optimalControlInstance_Sim < optimalControlInverseKKT & optimalControlDirectFmincon
    properties
        % ioc/doc properties
         windowMethod = 'sliding';
         windowAdvance = 20; % window slides
         windowGap = 20; % inds between spline points
         windowLength = 61;
        
%        windowMethod = 'segments';
%        windowAdvance = 200; % window slides
%        windowGap = 20; % inds between spline points
%        windowLength = 201;
        
        % doc fmincon properties
        fminconDisplaySetting = 'iter-detailed'; % 'final' 'notify' 'off' 'iter-detailed'
        fminconMaxIter = 1000;
        fminconTolFun = deg2rad(1e-6);
        fminconTolX = deg2rad(1e-6);
        fminconTolCon = deg2rad(1e-6);
        
        pivotMethod = 'resnorm';
        
        % end effector ind
%         endEffInd = 3;
        
%         normalizeArray = [1 1 1];
%         normalizeArray = 1./ [      435.481959543731          80.1772064356161          220044.871439183];
 
        normalizeArray = [0.0209193477500000, 0.391311233500000, 0.000555224862500000];
    end
    
    methods
        function obj = modifyModel(obj)
            % update link length
            length = [0.4140    0.4388    0.0939];
            mass = [ 7.2288   18.5238   10.6926];
            
            com{1} = [ 0.1606;     0.0583;          0];
            com{2} = [ 0.1764;     0.0683;          0];
            com{3} = [ 0.0223;     0.0141;          0];
            
            
            inert{1} = [    1.2387         0         0
                0         0         0
                0         0    1.2387];
            
            inert{2} = [    3.5671         0         0
                0         0         0
                0         0    3.5671];
            
            inert{3} = [  0.0943         0         0
                0         0         0
                0         0    0.0943];
            
            obj.featureSet.model.transforms(4).t(1:3, 4) = [length(1) 0 0]';
            obj.featureSet.model.transforms(6).t(1:3, 4) = [length(2) 0 0]';
            obj.featureSet.model.transforms(8).t(1:3, 4) = [length(3) 0 0]';
            
            for i = 1:3
                obj.featureSet.model.bodies(i).com = com{i};
                obj.featureSet.model.bodies(i).I   = inert{i};
                obj.featureSet.model.bodies(i).m   = mass(i);
            end
        end
        
        function obj = setupCostFunctions(obj)
            % initialize the cost function structure
            emptyCell = cell(3, 1);
            obj.costFunctionStruct = struct(...
                'name', {emptyCell{:}},...
                'function', []);
            
            cost_fct_calc_sq = @(feature, normalize)  normalize * sum(sum(feature .^ 2)) / (size(feature, 1));

            obj.costFunctionStruct(1).name = 'ddq';
            obj.costFunctionStruct(1).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.ddq, normalizeSet.ddq);
            obj.costFunctionStruct(2).name = 'ddx';
            obj.costFunctionStruct(2).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.ddx, normalizeSet.ddx);
            obj.costFunctionStruct(3).name = 'tau';
            obj.costFunctionStruct(3).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.tau, normalizeSet.tau);
        end
        
        function obj = setupCostNormalization(obj)
            obj.normalizeCostSet = rlFeatureSet_ioc();
            
            arrayLength = 1;
            obj.normalizeCostSet.ddq = obj.normalizeArray(1)*ones(1, arrayLength);
            obj.normalizeCostSet.ddx = obj.normalizeArray(2)*ones(1, arrayLength);
            obj.normalizeCostSet.tau = obj.normalizeArray(3)*ones(1, arrayLength); 
        end
        
        function obj = setupTimeArray(obj)
            switch obj.setupMode
                case 'docGen'
                otherwise       
            end
            
             obj.splineDt = 0.01;
             obj.windowTime = (0:obj.windowLength-1)*0.01;
        end
        
        function obj = setupConstraintPoints(obj)
            switch obj.setupMode
                case 'docGen'
                    splineViaPoints.inds = [1 ...
                        ceil(1*numel(obj.windowTime)/2) ...
                        numel(obj.windowTime)];
                    
                    splineViaPoints.q = [  1.376368005804248   0.072353981324302   1.013259948798985
                        1.203448266374773   0.915917831047703  -0.120610705719433
                        1.376368005804248   0.072353981324302   1.013259948798985];
                    splineViaPoints.dq = zeros(3, 3);
                    splineViaPoints.ddq = zeros(3, 3);
                    
                otherwise
                    splineViaPoints.inds = [1 ...
                        41 ...
                        61];
                    
                    splineViaPoints.inds = [1 ...
                        ceil(1*numel(obj.windowTime)/2) ...
                        numel(obj.windowTime)];

                    splineViaPoints.q = obj.featureSet.q(splineViaPoints.inds, :);
                    splineViaPoints.dq = obj.featureSet.dq(splineViaPoints.inds, :);
                    splineViaPoints.ddq = obj.featureSet.ddq(splineViaPoints.inds, :);
            end
            splineViaPoints.time = obj.windowTime(splineViaPoints.inds);
            
            obj.constPoints = splineViaPoints;
        end
        
        function splineViaPoints = setupSplineViaPoints(obj)
            switch obj.setupMode
                case 'docGen' 
                    obj.windowLength = 201;
                    
                otherwise
                   
            end
            
            splineViaPoints.inds = 1:obj.windowGap:obj.windowLength;
            splineViaPoints.time = obj.windowTime(splineViaPoints.inds);
            
            splineViaPoints.q = obj.featureSet.q(splineViaPoints.inds, :);
            splineViaPoints.dq = obj.featureSet.dq(splineViaPoints.inds, :);
            splineViaPoints.ddq = obj.featureSet.ddq(splineViaPoints.inds, :);
            
%             obj.splinePoints = splineViaPoints;
        end
    end
end