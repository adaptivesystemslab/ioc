classdef NNClassifier < AClassifier
% QDA Classifier class
    properties(SetAccess=protected)
       net_params = [10 10 10]; %Default 3 hidden layers 10 neurons each
       net = [];
       classes = []; 
    end
    
    methods
        
        function obj = NNClassifier(net_params)
           obj.net_params = net_params;
        end
        
        function init(obj,varargin)
        % Initialization function for NN Classifier
            obj.net_params = varargin{1};
        end
        
        function error = train(obj,input,output)
            obj.net = feedforwardnet(obj.net_params);
            
            %Get unique classes so we change to this number of outputs
            obj.classes = unique(output);
            new_output = zeros(numel(output),numel(obj.classes));
            for i=1:numel(obj.classes)
                new_output(output==obj.classes(i),i) = 1;
            end
            
%             obj.net = train(obj.net,input',new_output');
            obj.net = trainrp(obj.net,input',new_output');
            
            out = obj.classify(input);
            error = sum(output ~= out)/numel(output);
            
        end
        
        function out = classify(obj,input)
            out_tmp = obj.net(input')';
            out = zeros(size(out_tmp,1),1);
            for i=1:size(out,1)
                out(i) = obj.classes((out_tmp(i,:)==max(out_tmp(i,:))));
            end
        end
        
        function obj_copy = copy(obj)
           obj_copy = NNClassifier(obj.net_params); 
        end
        
    end    
end