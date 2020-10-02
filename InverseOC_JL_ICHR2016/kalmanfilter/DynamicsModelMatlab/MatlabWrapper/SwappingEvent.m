classdef SwappingEvent < event.EventData
    properties
        %Sensor Handles that swap if indeed both are swapping
        sensor;
        %Measurements that caused the swap
        mes;
        %Measures from match matrix
        measure;
        %Sensor Index
        sensor_index;
        %Input measurement index
        mes_index;
        
        iter_index;
    end
    methods
        function obj = SwappingEvent(sensor,mes,sensor_index,mes_index,measure,iter_index)
            obj.sensor = sensor;
            obj.mes = mes;
            obj.sensor_index = sensor_index;
            obj.mes_index = mes_index;
            obj.measure = measure;
            obj.iter_index = iter_index;
        end
    end
end