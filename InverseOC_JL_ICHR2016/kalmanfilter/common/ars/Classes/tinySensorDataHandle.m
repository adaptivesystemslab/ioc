classdef tinySensorDataHandle < imuDataHandle
    
    methods
        function obj = tinySensorDataHandle(varargin)
            obj = obj@imuDataHandle(varargin{1},varargin{2},varargin{3});
        end
        
    end
    
    
    
    methods (Access=protected)
       function obj = loadHeader(obj)
          xDoc = xmlread(obj.filepathHeader);
          
          name_item = xDoc.getElementsByTagName('mMyName');
          obj.name = char(name_item.item(0).getFirstChild.getData);
          sampling_rate_item = xDoc.getElementsByTagName('mSamplingRate');
          obj.samplingRate = str2num(sampling_rate_item.item(0).getFirstChild.getData);
          obj.dt = 1/obj.samplingRate;
           
       end
        
        
    end
end