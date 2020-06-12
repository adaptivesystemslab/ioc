classdef MissingEvent < event.EventData
    properties
        %Sensor Handle for which data went missing
        sensor;
        %Closest measurement to the missing data
        mes;
        %Measure from match matrix
        measure;
        %Sensor Index
        sensor_index;
        %Input measurement index
        mes_index;
        
        iter_index;
    end
    methods
        function obj = MissingEvent(sensor,mes,sensor_index,mes_index,measure,iter_index)
            obj.sensor = sensor;
            obj.mes = mes;
            obj.sensor_index = sensor_index;
            obj.mes_index = mes_index;
            obj.measure = measure;
            obj.iter_index = iter_index;
        end
    end
end