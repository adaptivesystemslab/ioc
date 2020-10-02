classdef IsoMapTransform < ATransform
% Pricipal Component Analysis transform class

    properties (SetAccess=protected)
        %The different components
        mapping = [];
        %How many components are used (default 1)
        degree = 1;
    end

    methods    
        function obj = IsoMapTransform(degree)
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
        % Train PCA Transform
        
            %no error in PCA
            [~, obj.mapping] = compute_mapping(input,'Isomap',obj.degree);
            error = 0;
        end
        
        function out = apply(obj,input)
            out = out_of_sample(input,obj.mapping);
        end
    end
end