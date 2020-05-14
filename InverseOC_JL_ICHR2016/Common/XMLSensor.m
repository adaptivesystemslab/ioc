classdef XMLSensor < handle
    %Sensor class to load stuff from .header and .data sensor files
    
    properties (SetAccess = private)
       name = [];
       samplingRate = [];
       date_collected = [];
       
       data = [];
    end
    
   methods
       function obj = XMLSensor(path) 
       
           %Open the .header file
           xDoc = xmlread(path);
           
           xRoot = xDoc.getDocumentElement;
           
           obj.name = char(xRoot.getAttribute('name'));
           obj.samplingRate = str2double(char(xRoot.getAttribute('samplingRate')));
           obj.date_collected = datenum(char(xRoot.getAttribute('date')),'yyyy_mm_dd_HH_MM_SS');
           
           obj.data = parseCSV(strrep(path,'.header','.data'));
           
       end
       
       
       function rotation = rotate2Grav(obj,grav_axis,window_size,modify)
            %Finds a window of points with smallest variance and uses it to
            %calculate a rotation such that T*data = norm(data) in grav_axis. 
            %Data is 3xN, grav axis is 'X' 'Y' 'Z' '-X' '-Y' '-Z'
            %If MODIFY is set to 1 it will rotate the acceleration and gyro
            %data for this sensor
            
            dat = [obj.data.AccelerometerXCalibrated';...
                        obj.data.AccelerometerYCalibrated';...
                        obj.data.AccelerometerZCalibrated'];
            
                switch grav_axis
                    case 'X'
                        grav = [1 0 0];
                    case '-X'
                        grav = [-1 0 0];
                    case 'Y'
                        grav = [0 1 0];
                    case '-Y'
                        grav = [0 -1 0];
                    case 'Z'
                        grav = [0 0 1];
                    case '-Z'
                        grav = [0 0 -1];
                    otherwise
                        grav = zeros(1,3);
                        axis = find(abs(mean(dat,2))==max(abs(mean(dat,2))),1);
                        grav(axis) = sign(dat(axis,1));
                end
                

                
                best_var = inf;
                best_index = 1;
                N = 1000;
                if N > (size(dat,2)-window_size)
                   N = size(dat,2)-window_size;
                end
                for i=1:N

                    variance = mean([var(dat(1,i:i+window_size))...
                        var(dat(2,i:i+window_size))...
                        var(dat(3,i:i+window_size))]);

                    if(variance < best_var)
                        best_var = variance;
                        best_index = i;
                    end
                end

                m = mean(dat(:,best_index:best_index+window_size),2);
                r0 = vrrotvec(grav*norm(m),m);
                rotation = vrrotvec2mat(r0)';
                
                
                
                if exist('modify','var') 
                   if modify
                       dat = rotation*dat;
                       obj.data.AccelerometerXCalibrated = dat(1,:)';
                       obj.data.AccelerometerYCalibrated = dat(2,:)';
                       obj.data.AccelerometerZCalibrated = dat(3,:)';
                       
                       dat = [obj.data.GyroscopeXCalibrated';...
                        obj.data.GyroscopeYCalibrated';...
                        obj.data.GyroscopeZCalibrated'];
                       
                       dat = rotation*dat;
                       
                       obj.data.GyroscopeXCalibrated = dat(1,:)';
                       obj.data.GyroscopeYCalibrated = dat(2,:)';
                       obj.data.GyroscopeZCalibrated = dat(3,:)';
                   end
                end
            end
       
            function offset = gyroOffset(obj,window_size,modify)
                %Calculates the gyro offset, if modify is 1 it will
                %subtract it from the data.
                
                dat = [obj.data.GyroscopeXCalibrated';...
                        obj.data.GyroscopeYCalibrated';...
                        obj.data.GyroscopeZCalibrated'];
                
                best_var = inf;
                best_index = 1;
                N = 1000;
                if N > (size(dat,2)-window_size)
                   N = size(dat,2)-window_size;
                end
                for i=1:N

                    variance = mean([var(dat(1,i:i+window_size))...
                        var(dat(2,i:i+window_size))...
                        var(dat(3,i:i+window_size))]);

                    if(variance < best_var)
                        best_var = variance;
                        best_index = i;
                    end
                end

                offset = mean(dat(:,best_index:best_index+window_size),2);
                
                if exist('modify','var')
                   if modify
                      dat = dat - repmat(offset,1,numel(dat(1,:))); 
                       obj.data.GyroscopeXCalibrated = dat(1,:)';
                       obj.data.GyroscopeYCalibrated = dat(2,:)';
                       obj.data.GyroscopeZCalibrated = dat(3,:)';
                   end
                end
                
            end
       
   end
    
    
    
end