classdef SLClassifier < AClassifier
   %This is the Support Vector Machine classifier
   %Uses matlab built in SVM functionality
   
    properties(Access=private)
       SVMmodel = [];
       modelParam = [];
       
       downsampleUpperBound = 1500; % the amount of points used for parameter searching
       applyZeroMean = 0; % 1 = remove mean and apply unity covar before training and testing
       meanOffset = [];
       stdDevOffset = [];
       
       %Default Kernel is RBF
       meanPhase = [];
       maxPhase = [];
       minPhase = [];
       
       ignoreVal = -Inf;
    end
    
    methods
        function init(obj,varargin)
            % set up kernel and tuning variables for the SVM training
            
        end
        
        function input = preprocessData(obj, input)
            
        end
        
        function error = train(obj,input, output)
            % this function performs a simple parameter search, then trains
            % the SVM based on the optimal parameter
            input = preprocessData(obj, input);
            
            % determine the range of the 'segment' points
            % first dof is the phase, second dof is the radius
            % pull out all points that is in the right quad, then find the
            % bounds            
            [meanPhase(1), maxPhase(1), minPhase(1)] = obj.meanRange(input(:, 1), input(:, 2), output,  0,    pi/2);
            [meanPhase(2), maxPhase(2), minPhase(2)] = obj.meanRange(input(:, 1), input(:, 2), output,  pi/2, pi);
            [meanPhase(3), maxPhase(3), minPhase(3)] = obj.meanRange(input(:, 1), input(:, 2), output, -pi,  -pi/2);
            [meanPhase(4), maxPhase(4), minPhase(4)] = obj.meanRange(input(:, 1), input(:, 2), output, -pi/2, 0);
            
            % rounding it off
            maxPhase(2) = pi;
            maxPhase(4) = 0;
            minPhase(1) = 0;
            minPhase(3) = -pi;
            
            obj.meanPhase = meanPhase;
            obj.maxPhase = maxPhase;
            obj.minPhase = minPhase;
        end
        
        function out = classify(obj,input)
            out = zeros(size(input, 1), 1);
            for ind = 1:size(input, 1)
                currPhase = input(ind, 1);
                
                lowerBound = currPhase > obj.minPhase;
                upperBound = currPhase < obj.maxPhase;
                
                if sum(and(lowerBound, upperBound))
                    out(ind) = 1;
                end
            end
        end

        function obj_copy = copy(obj)
           obj_copy = SVMClassifier();
           obj_copy.init; % linear would = 0 if the +1 didn't exist
        end
    end    
    
    methods(Static)
        function [meanPhase, maxPhase, minPhase] = meanRange(phaseVal, radVal, outputVal, minRange, maxRange)
            % find values that are within the range
            lowerBound = find(phaseVal >= minRange);
            upperBound = find(phaseVal <= maxRange);
            onBound = intersect(lowerBound, upperBound);
            
            phaseBound = phaseVal(onBound);
            outputBound = outputVal(onBound);
            
            phaseOfInterest = phaseBound(outputBound == 1);
            meanPhase = mean(phaseOfInterest);
            varPhase = var(phaseOfInterest);
            
            maxPhase = meanPhase + varPhase;
            minPhase = meanPhase - varPhase; 
%             maxPhase = max(phaseBound(outputBound == 1));
%             minPhase = min(phaseBound(outputBound == 1));
        end
    end
end