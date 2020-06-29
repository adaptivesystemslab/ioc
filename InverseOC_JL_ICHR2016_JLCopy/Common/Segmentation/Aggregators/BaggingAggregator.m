classdef BaggingAggregator < AAggregator
   
    properties(SetAccess=protected)
        
        aggregate_function = 1;
        
        classifier_set = [];
        
        %Unique classes
        classes = [];
        
        %Number of training sets
        n_sets = 10;
        
        %Number of samples per set, -1 means N/2
        n_samples = -1;
        
        %weights for LSE aggregation
        w = [];
        
        %Aggregating classifier
        aggregating_classifier = [];
    end
    
    
    methods
        
        function obj = BaggingAggregator(classifier)
            
            sc = superclasses(classifier);
            if ~find(ismember(sc,'AClassifier'))
                error('Input to BaggingAggregator Constructor must be a classifier');
            end
            
            obj.base_classifier = classifier;
        end
        
        function init(obj,varargin)
            %Initialization Input: 
            % Number of classifiers used (default 10)
            % Number of samples per classifier (default 2*n/3)
            % Aggregator Type (dafault majority voting)
                % type 1 : majority voting
                % type 2 : LSE-based weighting
                % type 2 : Classifier, in this case arg in 4 must be a
                % classifier
            
                if length(varargin) == 1
                    obj.n_sets = varargin{1};
                elseif length(varargin) == 3
                    obj.n_sets = varargin{1};
                    obj.n_samples = varargin{2};
                    obj.aggregate_function = varargin{3};
                end
            
            
            if nargin == 5
                obj.aggregating_classifier = varargin{4};
            end
        end
        
        function error = train(obj,input,output)
            
            %Use 2*n/3 sampling by default
            N = size(input,1);
            samples = obj.n_samples;
            if obj.n_samples == -1
                samples = floor(2*N/3);
            end
            
            %remember unique classes
            obj.classes = unique(output);
            
            %Create the calssifier set and train each classifier
            obj.classifier_set = cell(obj.n_sets,1);
            indexes = randi(N,obj.n_sets,samples);
            for i=1:obj.n_sets
                obj.classifier_set{i} = obj.base_classifier.copy();
                obj.classifier_set{i}.train(input(indexes(i,:),:),output(indexes(i,:,:)));
            end
            
            switch obj.aggregate_function
                case 1
                    %Majority based voting, dont need to do anything here
                case 2
                    %LSE-weights calculate A
                    A = zeros(N,obj.n_sets);
                    for j=1:obj.n_sets
                       A(:,j) = obj.classifier_set{j}.classify(input);
                    end
                    obj.w = A\output;
                case 3
                    %Train aggrigating classifier
                    agg_input = zeros(N,obj.n_sets);
                    for j=1:obj.n_sets
                       agg_input(:,j) = obj.classifier_set{j}.classify(input);
                    end
                    obj.aggregating_classifier.train(agg_input,output);                    
            end
            
            error = sum(obj.classify(input) ~= output)/N;
        end
        
        function out = classify(obj,input)
            
            N = size(input,1);
            A = zeros(N,obj.n_sets);
            out = zeros(N,1);
            for j=1:obj.n_sets
               A(:,j) = obj.classifier_set{j}.classify(input);
            end

            switch obj.aggregate_function
                case 1
                    %majority voting
                    for i=1:N
                       out(i,:) = mode(A(i,:)); 
                    end
                case 2
                    %LSE-weights only works for 2 class problem
                    for i=1:N
                        out(i,:) = interp1(obj.classes,obj.classes,obj.w'*(A(i,:)'),'nearest');
                    end
                    
                case 3
                    out = obj.aggregating_classifier.classify(A);
            end
        end
        
    end
end