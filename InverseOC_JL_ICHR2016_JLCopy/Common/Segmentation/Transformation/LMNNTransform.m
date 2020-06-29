classdef LMNNTransform < ATransform
% LMNN - uses the van der Maaten's DR MATLAB toolbox
% http://lvdmaaten.github.io/drtoolbox/

    properties (SetAccess=protected)
        %The different components
        mapping = [];
        %How many components are used (default 1)
        degree = 1;
        
        pcaModel = [];
        
         downsampleUpperBound = 2000; % so n p1 points, and n p0 points
    end

    methods    
        function obj = LMNNTransform(degree)
        % Constructor for PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = degree;
        end
        
        function init(obj,varargin)
        % Initialization of PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = varargin(1);
        end
        
        function error = train(obj,input,output)
            % Train PCA Transform to reduce error
            
              % apply downsample and save the training data
            if size(input, 1)/2 > obj.downsampleUpperBound
                % randomly select some points for downsampling
                p1points = find(output == 1); % don't want to assume the array is sorted
                p0points = find(output == 0);
                
                selectInd = randperm(length(p1points), obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(length(p0points), obj.downsampleUpperBound);
                
                X_tr = input([p1points(selectInd) p0points(selectInd2)], :);
                Y_tr = output([p1points(selectInd) p0points(selectInd2)], :);
            else
                X_tr = input;
                Y_tr = output;
            end
            
            obj.pcaModel = PCATransform(-1, 0.95);
            obj.pcaModel.train(input,output);
            inputRot = obj.pcaModel.apply(input);
        
            %no error in PCA
            [~, obj.mapping] = compute_mapping([Y_tr X_tr],'LMNN',obj.pcaModel.elbow);
            error = 0;
        end
        
        function out = apply(obj,input)
            inputRot = obj.pcaModel.apply(input);
            
            out = out_of_sample(inputRot,obj.mapping);
        end
    end
end