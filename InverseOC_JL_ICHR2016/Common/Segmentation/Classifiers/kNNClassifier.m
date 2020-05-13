classdef kNNClassifier < AClassifier
    % k-NN Classifier class
    properties(SetAccess=protected)
        trainingData = [];
        trainingLabel = [];
        
        mdl = [];
        
        k = 5;
        distanceMetric = 'euclidean';        
    end
    
    methods
        function obj = kNNClassifier(k, distanceMetric)
            % constructor 
            obj.init(k, distanceMetric);
        end
        
        function init(obj,k, distanceMetric)
            % Initialization function for QDA Classifier
            obj.k = k;
            obj.distanceMetric = distanceMetric;
        end
        
        function error = train(obj,input,output)
            obj.mdl = ClassificationKNN.fit(input, output, ...
                'NumNeighbors', obj.k, ...
                'Distance', obj.distanceMetric);
            
            error = [];
        end
        
        function [out, prob] = classify(obj,input,output,extra)
            
            out = zeros(size(input, 1), 1);
            for i = 1:size(input, 1)
                out(i) = obj.mdl.predict(input(i, :));
            end
            
            prob = ones(size(out));
        end
        
        function obj_copy = copy(obj)
           obj_copy = kNNClassifier(obj.k,obj.distanceMetric); 
        end
    end
end