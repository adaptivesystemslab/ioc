classdef EventHandler < handle
    
    properties
        miss_events = 0;
        swap_events = 0;
        mach_events = 0;
        
    end
    
    methods
        function matched_marker(obj,src,evtdata)
            obj.mach_events = obj.mach_events+1;
            
            disp(['Frame ' num2str(evtdata.iter_index) ': Marker ' num2str(evtdata.mes_index) ' Matched to ' evtdata.sensor.name]);
        end
        
        
        function swap_marker(obj,src,evtdata)
            obj.swap_events = obj.swap_events+1;
            
            if numel(evtdata.sensor) ==2
                disp(['Frame ' num2str(evtdata.iter_index) ': Marker ' evtdata.sensor(1).name ' Swapped With ' evtdata.sensor(2).name]);
            else
                disp(['Frame ' num2str(evtdata.iter_index) ': Marker ' evtdata.sensor(1).name ' Swapped ' num2str(evtdata.mes_index)]);
            end
        end
        
        function missing_marker(obj,src,evtdata)
            obj.miss_events = obj.miss_events+1;
            disp(['Frame ' num2str(evtdata.iter_index) ': Marker ' evtdata.sensor.name ' Missing']);
        end
    end
end