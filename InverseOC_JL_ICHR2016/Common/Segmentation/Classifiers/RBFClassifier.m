classdef RBFClassifier < AClassifier
% QDA Classifier class
    properties(SetAccess=protected)
       W = [];
       centers = [];
       cov = [];
       N = [];
       classes = [];
    end
    
    methods
        
        function obj = RBFClassifier(N)
           obj.N = N; 
        end
        
        function init(obj,varargin)
        % Initialization function for RBF Classifier
        % Input: Number of Radial Basis Functions
            obj.N = varargin{1};
        end
        
        function error = train(obj,input,output)
            opt = [0 0 1e-3 0 0 0 0 0 0 0 0 0 0 1000];
            mix = gmm(size(input,2), obj.N, 'spherical');
            mix = gmminit(mix,input,opt);
            [mix, ~, ~] = gmmem(mix,input,opt);

            obj.centers = mix.centres;
            obj.cov = mix.covars;

            G = exp(-dist(input,obj.centers').^2./ repmat((2.*obj.cov),size(input,1),1));
            
            %Get unique classes so we change to this number of outputs
            obj.classes = unique(output);
            new_output = zeros(numel(output),numel(obj.classes));
            for i=1:numel(obj.classes)
                new_output(output==obj.classes(i),i) = 1;
            end
            
            
            obj.W = linsolve([G ones(size(input,1),1)],new_output);
            out = obj.classify(input);
            error = sum(output ~= out)/numel(output);
        end
        
        function out = classify(obj,input)
                G = exp(-dist(input,obj.centers').^2./ repmat((2.*obj.cov),size(input,1),1));
                out_tmp = [G ones(size(input,1),1)]*obj.W;
                out = zeros(size(out_tmp,1),1);
                for i=1:size(out,1)
                    out(i) = obj.classes(out_tmp(i,:)==max(out_tmp(i,:)));
                end
        end
        
        function obj_copy = copy(obj)
           obj_copy = RBFClassifier(obj.N); 
        end
        
    end    
end