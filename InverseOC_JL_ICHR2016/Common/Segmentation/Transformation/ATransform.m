classdef ATransform < handle    
% Abstract Class for a transofrm 
% Subclass must implement init, train and apply
% the operator * (Obj*Input) can also be used to apply the transform
    
    methods(Abstract)
        init(obj,varargin)
        error = train(obj,input,output)
        out = apply(obj,input)        
    end
    methods(Sealed)
       function out = mtimes(obj,input)
       % Override the * operator so we can do Output = Transform*Input 
           
            out = apply(obj,input);
        end 
    end
end