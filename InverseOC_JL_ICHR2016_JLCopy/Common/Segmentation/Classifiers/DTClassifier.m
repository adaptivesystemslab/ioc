classdef DTClassifier < AClassifier
% FSM
    properties(SetAccess=protected)
       training = [];
       group = [];
       
       model = [];
       
       labelMapping = {};
    end
    
    methods
        function init(obj,varargin)
        % Initialization function for FSM

        end
        
        function error = train(obj,trainingData, trainingLabel)
          
            obj.model = ClassificationTree.fit(trainingData, trainingLabel);
            
            classifyLabel = obj.classify(trainingData);
            error = 1 - sum(trainingLabel == classifyLabel) / length(classifyLabel);
        end
        
        function [out, prob] = classify(obj,input)
            [out, prob] = predict(obj.model,input);
        end
        
        function obj_copy = copy(obj)
            obj_copy = finiteStateMachine();
        end
    end    
end