classdef LassoTransform < ATransform
    % Lasso dim selection code
    
    properties(SetAccess=private)
        
       lassoFeatures = [];
       degree = 2;
       
       meanOffset = [];
       stdDevOffset = [];
       
       applyZeroMean = 1; % 1 = remove mean and apply unity covar before training and testing
       downsampleUpperBound = 5000; % so n p1 points, and n p0 points
    end
        
    methods
        function obj = LassoTransform(varargin)
            obj.degree = varargin{1};
        end
        
        function init(obj,varargin)
            
        end
        
        function error = train(obj,X,Y)
            if size(X, 1) > obj.downsampleUpperBound
                % randomly select some points for downsampling
                halfXLength = floor(size(X, 1)/2);
                selectInd = randperm(halfXLength, obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(halfXLength, obj.downsampleUpperBound) + halfXLength;
                X = X([selectInd selectInd2], :);
                Y = Y([selectInd selectInd2], :);
            else

            end
            
            % The data shold be normalized before inserting into the 
            % system: zero mean and unity covariance. store the offset
            % covar to apply to the testing data
            if obj.applyZeroMean
                [X, obj.meanOffset, obj.stdDevOffset] = zeroMeanCovarSignal(X, 1, 1);
            else 
                
            end
            
            [B,FitInfo] = lasso(X,Y);
            B_abs = ceil(abs(B));
            B_sum = sum(B_abs > 0, 1);
            [closeVal, closeInd] = findClosestValue(obj.degree, B_sum, 'smaller'); % pulls out a value, but the array gets sorted
            BVectorSelected = find(B_sum == closeVal);
            
            obj.lassoFeatures = B_abs(:, BVectorSelected(1)) == 1;
        end
        
        function out = apply(obj,input)
            if obj.applyZeroMean
                % shift the data to 0 mean, 1 covar
                input = zeroMeanCovarSignal(input, 1, 1, obj.meanOffset, obj.stdDevOffset);
            else 
                
            end
            
            out = input(:, obj.lassoFeatures);
        end
    end
end