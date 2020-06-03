classdef QDAClassifier < AClassifier
% QDA Classifier class
    properties(SetAccess=protected)
       training = [];
       group = [];
        
    end
    
    methods
        function init(obj,varargin)
        % Initialization function for QDA Classifier
            
        end
        
        function error = train(obj,input,output)
            obj.training = input;
            obj.group = output;
            [class,error] = classify(input,obj.training,obj.group,'quadratic');
        end
        
        function out = classify(obj,input)
            out = classify(input,obj.training,obj.group,'quadratic');
        end
        
        function obj_copy = copy(obj)
            obj_copy = QDAClassifier();
        end
        
    end    
end