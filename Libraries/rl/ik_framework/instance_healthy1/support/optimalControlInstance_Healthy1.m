classdef optimalControlInstance_Healthy1 < optimalControlInverseKKT & optimalControlDirectFmincon
    properties
        % ioc/doc properties
         windowMethod = 'sliding';
         windowAdvance = 20; % window slides
         windowGap = 20; % inds between spline points
         windowLength = 61;
        
        
        % doc fmincon properties
        fminconDisplaySetting = 'iter-detailed'; % 'final' 'notify' 'off' 'iter-detailed'
        fminconMaxIter = 1000;
        fminconTolFun = deg2rad(1e-6);
        fminconTolX = deg2rad(1e-6);
        fminconTolCon = deg2rad(1e-6);
        
        pivotMethod = 'resnorm';
        
        normalizeArray = [1 1 1 1 1 1 1 1];
%         normalizeArray = 1./ [      435.481959543731          80.1772064356161          220044.871439183];
    end
    
    methods        
        function obj = setupCostFunctions(obj)
            % initialize the cost function structure
%             emptyCell = cell(8, 1);
%             obj.costFunctionStruct = struct(...
%                 'name', {emptyCell{:}},...
%                 'function', []);
            
            cost_fct_calc_sq  = @(feature, normalize)  normalize * sum(sum(feature .^ 2)) / (size(feature, 1)*size(feature, 2));
            cost_fct_calc_abs = @(feature, normalize)  normalize * sum(sum(abs(feature))) / (size(feature, 1)*size(feature, 2));

            obj.costFunctionStruct(1).name = 'ddq';
            obj.costFunctionStruct(1).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.ddq, normalizeSet.ddq);
            obj.costFunctionStruct(2).name = 'dddq';
            obj.costFunctionStruct(2).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.dddq, normalizeSet.dddq);
            obj.costFunctionStruct(3).name = 'dddx';
            obj.costFunctionStruct(3).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.dddx, normalizeSet.dddx);
            obj.costFunctionStruct(4).name = 'tau';
            obj.costFunctionStruct(4).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.tau, normalizeSet.tau);
            obj.costFunctionStruct(5).name = 'dtau';
            obj.costFunctionStruct(5).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.dtau, normalizeSet.dtau);
            obj.costFunctionStruct(6).name = 'ddtau'; % torque acceleration (Berret2011: 'effort')
            obj.costFunctionStruct(6).function = @(featureSet, normalizeSet) cost_fct_calc_sq(featureSet.ddtau, normalizeSet.ddtau);
            obj.costFunctionStruct(7).name = 'power'; % angular power (Berret2011: 'energy')
            obj.costFunctionStruct(7).function = @(featureSet, normalizeSet) cost_fct_calc_abs(featureSet.power, normalizeSet.power);
            obj.costFunctionStruct(8).name = 'ek';  % kinetic energy (modified from Berret2011's 'geodesic', which where ek = geodesic^2)
            obj.costFunctionStruct(8).function = @(featureSet, normalizeSet) cost_fct_calc_abs(featureSet.ek, normalizeSet.ek);
            
        end
        
        function obj = setupCostNormalization(obj)
            obj.normalizeCostSet = rlFeatureSet_ioc();
            
            arrayLength = 1;
            obj.normalizeCostSet.ddq = obj.normalizeArray(1)*ones(1, arrayLength);
            obj.normalizeCostSet.dddq = obj.normalizeArray(2)*ones(1, arrayLength);
            obj.normalizeCostSet.dddx = obj.normalizeArray(3)*ones(1, arrayLength);
            obj.normalizeCostSet.tau = obj.normalizeArray(4)*ones(1, arrayLength);
            obj.normalizeCostSet.dtau = obj.normalizeArray(5)*ones(1, arrayLength);
            obj.normalizeCostSet.ddtau = obj.normalizeArray(6)*ones(1, arrayLength);
            obj.normalizeCostSet.power = obj.normalizeArray(7)*ones(1, arrayLength);
            obj.normalizeCostSet.ek = obj.normalizeArray(8)*ones(1, arrayLength);
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