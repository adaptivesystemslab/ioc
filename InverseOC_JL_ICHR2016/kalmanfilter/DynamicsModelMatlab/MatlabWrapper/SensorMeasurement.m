classdef SensorMeasurement < handle
    properties
        type = 0;
        size = 0;
        mes = [];
    end
    
    methods
        function obj = SensorMeasurement(varargin)
            %Constructor that extracts measurement and type from given array
            %of sensors
            if nargin ~= 0
                %If our input is an array of sensor objects use their
                %measurements
                if(isa(varargin{1},'SensorCore'))
                    sensors = varargin{1};
                    n = size(sensors,1);
                    m = size(sensors,2);
                    obj(n,m) = SensorMeasurement;
                    for i=1:numel(sensors)
                        obj(i).type = sensors(i).binary_type;
                        obj(i).mes = sensors(i).measurement;
                        obj(i).size = numel(obj(i).mes);
                    end
                else
                    %Just m by n measurement array
                    obj(varargin{1},varargin{2}) = SensorMeasurement;
                end
                
            end
        end
        
        function value = getMesArray(obj)
           m = size(obj,1);
           n = numel(vertcat(obj(1,:).mes));
           value=zeros(m,n);
           for i=1:m
               value(i,:) = vertcat(obj(i,:).mes);
           end
        end
        
        
        function setMesArray(obj,value)
            %Assign Measurement to array of SensorMeasurement Objects
            
            %Single SensorMeasurement assignment
            if numel(obj)==1
               obj.mes = value;
               obj.size = numel(value);
               return;
            end
            
            %Array Assignment
            [m, n] = size(obj);
            for i=1:m
                mes_col_indx = 1;
                for j=1:n
                    obj(i,j).mes = value(i,mes_col_indx:mes_col_indx+obj(i,j).size-1)';
                    mes_col_indx = mes_col_indx+obj(i,j).size;
                end
            end
        end
    end
end