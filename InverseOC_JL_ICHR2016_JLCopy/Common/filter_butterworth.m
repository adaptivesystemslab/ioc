classdef filter_butterworth
    % PURPOSE - provides a Butterworth filter
    % Accepts 1-dim or 2-dim (vertical vectors) data
    
    % type          filter type: low, high, band
    % pass          how many passes: single, double
    % cutoff        cutoff frequency
    % sampleRate    sample rate of the data
    % order         order of Butterworth.
    % zi            initial values for the filter
    % padLength     if non-zero, denotes the amount of pre-padding that
    %               should be added to the data to reduce rise time issues
    
    % output        filtered data
    % zf            final condition of the butterworth filter
    
    properties(SetAccess=protected)
        % filter parameters
        type = 'low';
        pass = 'single';        
        filterFreq = 0.04;
        order = 2;        

        % internal filter values
        a = [];
        b = [];
    end
    
    properties(SetAccess=public)
        padLength = 0;
    end
    
    methods        
        function obj = filter_butterworth(type, pass, filterFreq, order, padLength)
            % set up the filter parameters
            if exist('type', 'var') || ~empty(type)
                obj.type = type;
            end
            
            if exist('pass', 'var') || ~empty(pass)
                obj.pass = pass;
            end
            
            if exist('filterFreq', 'var') || ~empty(filterFreq)
                obj.filterFreq = filterFreq;
            end
            
            if exist('order', 'var') || ~empty(order)
                obj.order = order;
            end
            
            if exist('padLength', 'var') || ~empty(padLength)
                obj.padLength = padLength;
            end
            
             % generate the filter parameters
             [obj.b, obj.a] = butter(obj.order, obj.filterFreq, obj.type);
        end
        
        function [output, zf] = apply(obj, filterVal, zi)
            dataDof = size(filterVal, 2);
                
            if ~exist('zi', 'var') || isempty(zi)
                zi = zeros(obj.order, dataDof);
            end
            
            % generate temp variable
            output = zeros(size(filterVal));
            zf = zeros(obj.order, dataDof);
            
            for ind_dof = 1:dataDof
                arrayToFilter = filterVal(:, ind_dof);
                
                if obj.padLength
                    % insert a padding to 'wind up' the filter
                    arrayToFilter = padarray(arrayToFilter, obj.padLength, arrayToFilter(1), 'pre');
                end
                
                [dofOutput, dofZf] = obj.filterFunction(arrayToFilter, zi(:, ind_dof));
                
                if obj.padLength
                    % then crop out the 'wind up' rise time
                    dofOutput = dofOutput(obj.padLength+1:end);
                end
                
                output(:, ind_dof) = dofOutput;
                zf(:, ind_dof) = dofZf;
            end
        end
        
        function [output, zf]  = filterFunction(obj, filterVal, zi)
            % apply the filter
            
            if strcmpi(obj.pass, 'singlePass')
                [output, zf] = filter(obj.b, obj.a, filterVal, zi);
            elseif strcmpi(obj.pass, 'doublePass')
                [output, zf]  = filtfilt(obj.b, obj.a, filterVal, zi);
            elseif strcmpi(obj.pass, 'noFilter')
                output = filterVal;
                zf = zi;
            end
        end
    end
end