classdef ADataHandle < matlab.mixin.Copyable
% ADatahandle abstract data handle for the individual data saving and 
% loading capabilities to inherent. This class will hold the capability to:
% - load the data given a file path
% - store the data in a structure that makes sense for the data on hand
% - provide a plotting capability to visualize the data 
%
% the contents of this class will be used in conjunction with
% exerciseDataset, which will contain an instance of this class, so
% demographic information such as subject number and exercise type should
% already be available there
%
% if the nature of the data is such that there are multiple instances of a
% given data, such as the IMU devices, use this class as an array
%
% last updated: March 25, 2014

    properties (SetAccess = protected)
        filepathHeader % the filepath to the header
        filepathData   % the filepath to the data itself
    end
    
    properties (SetAccess = public)
        % the individual data structure will depend on the underlying
        % data itself
    end
    
    methods (Abstract) % user will access these data
        load(obj,varargin) % load header and datafrom file
        write(obj,varargin) % write header and data to file
        h = plot(obj,varargin) % plot data
    end
    
    methods (Access = private) % the subfuctions should be hidden from the user
        loadHeader(obj,varargin)  % load the contents of the header from file
        loadData(obj,varargin)    % load the contents of the data from file
        writeHeader(obj,varargin) % write the contents of the header to file
        writeData(obj,varargin)   % write the contents of the data to file
    end
end