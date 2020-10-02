classdef AClassifier < handle
% Abstract Classifier class 
% Subclasses must implement init, train, classify
    
    methods(Abstract)
        init(obj,varargin)
        error = train(obj,input,output)
        out = classify(obj,input) 
        obj_copy = copy(obj)
    end    
end