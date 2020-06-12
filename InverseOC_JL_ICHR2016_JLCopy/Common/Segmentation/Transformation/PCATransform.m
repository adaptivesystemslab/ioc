classdef PCATransform < ATransform
% Pricipal Component Analysis transform class

    properties 
        %The different components
        U = [];
        m = [];
        %How many components are used (-1 means calculate the elbow of eigenvalues)
        degree = -1;
        elbow = [];
        
        threshold = 0;
    end

    methods    
        function obj = PCATransform(degree, threshold)
        % Constructor for PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = degree;
            
            if ~exist('threshold', 'var')
                obj.threshold = 0.8;
            else
                obj.threshold = threshold;
            end
        end
        
        function init(obj,varargin)
        % Initialization of PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = varargin(1);
        end
        
        function error = train(obj,input,output)
        % Train PCA Transform
        
            %no error in PCA
            error = 0;
            obj.m = mean(input);
            input = input-repmat(obj.m,numel(output),1);
            [U_ S ~]=svd(input','econ');
            obj.U = U_;
            
            %Calculate Elbow if we chose to
            if obj.degree == -1
                eig_vals = diag(S);
                eig_vals = eig_vals/sum(eig_vals);
                
                sum_eig = 0;
                index = 0;
                while sum_eig < obj.threshold
                    sum_eig = sum_eig + eig_vals(index+1);
                    index = index +1;
                end
                obj.elbow = index;
            end
        end
        
        function out = apply(obj,input)
%             input = bsxfun(@minus, input, obj.mapping.mean);
%             out = input * obj.mapping.M;
            if(obj.degree == -1)
                out = (obj.U(:,1:obj.elbow)'*(input-repmat(obj.m,size(input,1),1))')';
            else
                out = (obj.U(:,1:obj.degree)'*(input-repmat(obj.m,size(input,1),1))')';
            end
        end
    end
end