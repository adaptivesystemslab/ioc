classdef finiteStateMachine < AClassifier
% FSM
    properties(SetAccess=protected)
       training = [];
       group = [];
       
       fsm = [];

       %     auxState.rawRest = []; % 1 resting
       %     auxState.rawFlex = []; % 2 moving towards peak
       %     auxState.rawPeak = []; % 3 around peak
       %     auxState.rawExt = [];  % 4 coming back from peak
       labels = [1 0 1 0];
    end
    
    methods
        function init(obj,varargin)
        % Initialization function for FSM

        end
        
        function error = train(obj,input,output)
            
            if isstruct(input)
                trainingData = [input.rawRest; input.rawFlex; input.rawPeak; input.rawExt];

                trainingLabel = [obj.labels(1)*ones(size(input.rawRest, 1), 1); ...
                    obj.labels(2)*ones(size(input.rawFlex, 1), 1); ...
                    obj.labels(3)*ones(size(input.rawPeak, 1), 1); ...
                    obj.labels(4)*ones(size(input.rawExt, 1), 1)];
            else
                trainingData = input;
                trainingLabel = output;
            end
            
            if ~isempty(trainingData)
                obj.fsm = ClassificationTree.fit(trainingData, trainingLabel);
            end
        end
        
        function [out, prob] = classify(obj,input)
            [out, prob] = predict(obj.fsm,input);
        end
        
        function obj_copy = copy(obj)
            obj_copy = finiteStateMachine();
        end
    end    
end