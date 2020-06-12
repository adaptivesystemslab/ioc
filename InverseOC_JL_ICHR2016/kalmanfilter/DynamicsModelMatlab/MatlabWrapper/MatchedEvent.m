classdef MatchedEvent < event.EventData
    properties
        sensor;
        mes;
        measure;
        sensor_index;
        mes_index;
        iter_index;
    end
    methods
        function obj = MatchedEvent(sensor,mes,sensor_index,mes_index,measure,iter_index)
            obj.sensor = sensor;
            obj.mes = mes;
            obj.sensor_index = sensor_index;
            obj.mes_index = mes_index;
            obj.measure = measure;
            obj.iter_index = iter_index;
        end
    end
end