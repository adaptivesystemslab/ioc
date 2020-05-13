classdef AAggregator < handle
% Abstract Aggregator class 
% Subclasses must implement init, train, classify
    
    properties(Access = protected)
       base_classifier = []; 
    end

    methods(Abstract)
        init(obj,varargin)
        error = train(obj,input,output)
        out = classify(obj,input) 
    end    
end