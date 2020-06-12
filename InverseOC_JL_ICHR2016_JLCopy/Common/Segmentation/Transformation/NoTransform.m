classdef NoTransform < ATransform
    % No transformation 

    properties (SetAccess=protected)
        
    end

    methods    
        function obj = NoTransform()
        % Constructor for no transform
        end
        
        function init(obj,varargin)
        % Initialization of no transform
        
        end
        
        function error = train(obj,input,output)
       
        end
        
        function out = apply(obj,input)
%             input = bsxfun(@minus, input, obj.mapping.mean);
%             out = input * obj.mapping.M;
            
            out = input;
        end
    end
end