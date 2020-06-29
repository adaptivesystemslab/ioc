classdef EKF_Resize < EKF
    %Here we have EKF based on a kinematic model that is automatical resized
    %if sensors are missing
    properties(SetAccess = protected)
        %Flag if the sensor is missing or not (true = visible, false = missing)
        visible_indexes = logical.empty();
        visible_sensors = logical.empty();
        
        mes_sizes = [];
        mes_starts = [];
        sens_types = [];
        
        %This keeps track of the assignment indeces between incoming
        %measurements and sensors to trigger swap, missing, and found
        %callbacks.
        %For example we will have a structure with 
        %type = 1 : markers
        %mes = each row corresponds to marker attached to model 
        %       the entry is the index of previously used incoming
        %       measurement
        mes_assignment = SensorMeasurement.empty();
        
        %This is the final selected columns from the big K and H
        %Resizable EKF basically pulls these out in the update step
        selectedcols = [];
        
        %Threshold for detection of sensor gone missing
        %[type of sensor, threshold]
        sensor_missing_thresh = [1 100;2 1; 4 100; 17 100; 19 5; 7 1; 3 1;6 100; 8 100; 512 1000; 14 100];
        
        model_handle = [];        
    end
    
    events
        %Callbacks for events
        missing_callback;
        matched_callback;
        swap_callback;
    end
    
    methods
        
        function obj = EKF_Resize(model)
           %Constructor that sets everything visible 
           obj.model_handle = model;
           %Save the sizes of the measurements 
           obj.mes_sizes = cellfun(@numel,{obj.model_handle.sensors.measurement});
           %The start of each measurement
           obj.mes_starts = cumsum(obj.mes_sizes)-obj.mes_sizes+1;
           obj.sens_types = [model.sensors.binary_type]';
           %We dont know any sensor -> measurement assignments
           unique_types = unique(obj.sens_types);
           obj.mes_assignment = SensorMeasurement(numel(unique_types),1);
           for i=1:numel(unique_types)
               obj.mes_assignment(i).type = unique_types(i);
               obj.mes_assignment(i).mes = zeros(numel(...
                   obj.sens_types(obj.sens_types == unique_types(i))),1);
           end
        end
        
        function z = makeMeasure(obj,x)
            %Predict measurement from state
            dof = numel(obj.model_handle.joints);
            obj.model_handle.position = x(1:dof);
            obj.model_handle.velocity = x(dof+1:dof*2);
            obj.model_handle.acceleration = x(dof*2+1:dof*3);
            obj.model_handle.forwardPosition();
            obj.model_handle.forwardVelocity();
            obj.model_handle.forwardAcceleration();
            %Generate measurements
            z = SensorMeasurement(obj.model_handle.sensors);
        end
        
        %Run EKF Iteration
        function [cur_ass, cur_ass_out] = run_iteration(obj,u,z,matches)
            % Run Magical EKF Iteration 
            % Inputs :
            %       u: Delta Time
            %       z: SensorMeasurement object array
            % matches: Known correspondances between sensors and
            %          measurements for example
            %          matches = [1, 2] means that first attached sensor
            %          is assosiated with z(2) measurement            
            
            obj.iterationStep = obj.iterationStep + 1;
            
            obj.x_predict = obj.stateUpdate(obj.state,u);              %State Update Equation
            obj.A = obj.makeA(obj.x_predict,u);                        %State Update Jacobian
            obj.H = obj.makeH(obj.x_predict);                          %Measurement Jacobian
            
            obj.P_predict = obj.A*obj.covariance*obj.A'+obj.process_noise;%partial update
            obj.inov_covariance = obj.H*obj.P_predict*obj.H' + obj.MakeObservationNoise();              %cross covariance
            
            %This is used by the basic makeK method as well as for
            %measurement matching
            obj.inov_covariance_inv = inv(obj.inov_covariance);         %Compute inverse of the inovation covariance
            
            obj.z_predict = obj.makeMeasure(obj.x_predict);          %Measurement Prediction
            
            matchMatrices = obj.buildMatchMatrices(z,obj.z_predict,obj.inov_covariance_inv);
            
            %Here we do the matching and trigger callbacks for swapping,
            %missing, and matching sensors
            % Go through each sensor type
            cur_ass_out = SensorMeasurement(1, numel(matchMatrices));
            for i=1:numel(matchMatrices)
                %Force Known Matches, this can be better done so we dont 
                %spend time on already matched sensors
                %TODO, maybe figure out a way without sorting so faster
                if nargin == 4
%                     matches_type = matches(matches_subarray,:); % VJ orig
                    matches_subarray = find(obj.sens_types(matches(:,1)) == matchMatrices(i).type); % JL added 4 new lines
                    match_left = matches(matches_subarray, 1) - matches(matches_subarray(1), 1) + 1;
                    match_right = matches(matches_subarray, 2);
                    matches_type = [match_left match_right];
                    if(~isempty(matches_type))
                        %Get input mesurements of the same type
                        mes_indx = find([z.type] == matchMatrices(i).type);
                        %Find indexes
                        [C,ia,ib] = intersect(matches_type(:,2),mes_indx,'stable');
                        ia = matches_type(ia, 1); % JL added 
                        matchMatrices(i).mes(sub2ind(size(matchMatrices(i).mes),ia,ib)) = -inf;
                    end
                end
                mes = z([z.type] == matchMatrices(i).type);
                %These are the current assignments for this type of
                %sensor
                cur_ass = obj.mes_assignment([obj.mes_assignment.type] == matchMatrices(i).type);
                cur_ass_out(i) = cur_ass;
                %Go through each sensor attached to the model
                for j=1:size(matchMatrices(i).mes,1)
                    
                    %for this sensor the best assignment
                    [min_logP_sens, min_mark_indx] = min(matchMatrices(i).mes(j,:));
                    
                    %For the sens selected best assignment, the best sensor
                    [min_logP_mark, min_sens_indx] = min(matchMatrices(i).mes(:,min_mark_indx));
                    
                    threshold = obj.sensor_missing_thresh(obj.sensor_missing_thresh(:,1)==matchMatrices(i).type,2);
                    %We do not have an assignement for this sensor before but now
                    %we found a measurement above the threshold of logP
                    %and it is the lowest for all sensors
                    if(min_logP_sens < threshold && cur_ass.mes(j) == 0 &&...
                            min_logP_sens == min_logP_mark && ...
                            isempty(find(cur_ass.mes == min_mark_indx,1)))
                        
                        matched_event_obj = MatchedEvent(obj.model_handle.sensors(j),...
                            mes(min_mark_indx),j,min_mark_indx,min_logP_mark,obj.iterationStep);
                        notify(obj,'matched_callback',matched_event_obj);
                        %Save the match
                        cur_ass.mes(j) = min_mark_indx;
                        
                        %disp('MARKER MATCHED');
                    elseif(min_logP_sens < threshold && cur_ass.mes(j) ~= min_mark_indx && min_logP_sens == min_logP_mark)

                        %Check if the marker previously assigned to this
                        %sensor is the best assignment for the sensor that
                        %holds the current best marker assigned to it.
                        old_sensor_indx = find(cur_ass.mes == min_mark_indx,1);
                        if(~isempty(old_sensor_indx))
                            [min_logP_old_sens, min_new_mark_indx] = min(matchMatrices(i).mes(old_sensor_indx,:));
                            if(min_new_mark_indx == cur_ass.mes(j))
                                %Actual swap occured between 2 already
                                %assigned markers
                                cur_ass.mes(j) = min_mark_indx;
                                cur_ass.mes(old_sensor_indx) = min_new_mark_indx;
                                
                                mes = z([z.type] == matchMatrices(i).type);
                                swap_event_obj = SwappingEvent(...
                                    [obj.model_handle.sensors(j) obj.model_handle.sensors(old_sensor_indx)] ,...
                                    [mes(min_mark_indx) mes(min_new_mark_indx)],...
                                    [j old_sensor_indx],...
                                    [min_mark_indx min_new_mark_indx],...
                                    [min_logP_mark min_logP_old_sens],obj.iterationStep);
                                notify(obj,'swap_callback',swap_event_obj);
                                
                                %disp(['SWAP BETWEEN: ' num2str([j old_sensor_indx])]);
                            else
                                cur_ass.mes(j) = min_mark_indx;
                                cur_ass.mes(old_sensor_indx) = 0;
                                swap_event_obj = SwappingEvent(...
                                    [obj.model_handle.sensors(j) obj.model_handle.sensors(old_sensor_indx)] ,...
                                    [mes(min_mark_indx) mes(min_new_mark_indx)],...
                                    [j old_sensor_indx],...
                                    [min_mark_indx 0],...
                                    [min_logP_mark min_logP_old_sens],obj.iterationStep);
                                notify(obj,'swap_callback',swap_event_obj);
                            end
                        else
                            %The swap occured with a marker that has no
                            %assignment
                            swap_event_obj = SwappingEvent(...
                                obj.model_handle.sensors(j),...
                                mes(min_mark_indx),...
                                j,...
                                min_mark_indx,...
                                min_logP_mark,obj.iterationStep);
                            notify(obj,'swap_callback',swap_event_obj);
                            cur_ass.mes(j) = min_mark_indx;
                        end

                    elseif cur_ass.mes(j) ~=0 && (min_logP_sens > threshold ||...
                            min_logP_sens ~= min_logP_mark)
                        
                        mes = z([z.type] == matchMatrices(i).type);
                        missing_event_obj = MissingEvent(obj.model_handle.sensors(j),...
                            mes(min_mark_indx),j,min_mark_indx,min_logP_mark,obj.iterationStep);
                        notify(obj,'missing_callback',missing_event_obj);
                        cur_ass.mes(j) = 0;
                        
                        %disp('MARKER MISSING');
                    end
                end
            end
            
            %Convert everything correctly to matrices by picking correct
            %columns of K for matched assignment only
            obj.selectedcols = false(sum(obj.mes_sizes),1);
            DZ_vec = zeros(numel(vertcat(obj.z_predict.mes)),1);
            for i=1:numel(obj.mes_assignment)
               cur_type = obj.mes_assignment(i).type;
               
               cur_type_starts = obj.mes_starts(obj.sens_types == obj.mes_assignment(i).type);
               cur_type_starts = cur_type_starts(obj.mes_assignment(i).mes ~= 0);
               cur_type_sizes = obj.mes_sizes(obj.sens_types == obj.mes_assignment(i).type);
               cur_type_sizes = cur_type_sizes(obj.mes_assignment(i).mes ~= 0);
               
               z_pred = obj.z_predict([obj.z_predict.type] == cur_type);
               z_pred = z_pred(obj.mes_assignment(i).mes ~= 0);
               
               z_act = z([z.type] == cur_type);
               z_act = z_act(obj.mes_assignment(i).mes(obj.mes_assignment(i).mes ~= 0));
               
               for k=1:numel(cur_type_starts)
                    DZ_vec(cur_type_starts(k):cur_type_starts(k)+cur_type_sizes(k)-1) =...
                        z_act(k).mes - z_pred(k).mes;
                    obj.selectedcols(cur_type_starts(k):cur_type_starts(k)+cur_type_sizes(k)-1) = 1;
               end               
            end
            
            
            obj.dz = DZ_vec(obj.selectedcols);
            obj.H = obj.H(obj.selectedcols,:);
            obs_noise = obj.MakeObservationNoise();
            obs_noise = obs_noise(obj.selectedcols,obj.selectedcols);
            obj.inov_covariance =  obj.H*obj.P_predict*obj.H' + obs_noise;
            obj.K = obj.makeK();
            obj.state=obj.makeState();
            obj.covariance = obj.makeP(); %Updated Covariance Estimate
        end
        
        function logP = buildMatchMatrices(obj,z,z_predict,cov_inv)
            % Build the match matrices for each type of sensor
            % logP output will be SensorMeasurement array where each element
            % has a type and the matchmatrix corresponding to measurements
            % The match matrix is incoming measurements long and number of
            % sensors wide. 
            
            %Get unique type of measurement input
            types = unique([z.type]);
            logP = SensorMeasurement(numel(types),1);
            
            for type_indx = 1:numel(types)
                logP(type_indx).type = types(type_indx);
                %z_p = [z_predict([z_predict.type] == types(type_indx)).mes]';
                %z_a = [z([z.type] == types(type_indx)).mes]';
                %logP(type_indx).mes = inf(size(z_p,1),size(z_a,1));
                z_p = [z_predict([z_predict.type] == types(type_indx)).mes];
                z_a = [z([z.type] == types(type_indx)).mes];
                
                %We need to pull out appropriate covariance parts
                cs = obj.mes_starts(obj.sens_types == types(type_indx));
                ce = obj.mes_starts(obj.sens_types == types(type_indx)) + ...
                    obj.mes_sizes(obj.sens_types == types(type_indx))-1;
                
                %                 for i=1:size(z_a,1)
                %                     logP(type_indx).mes(:,i) = pdist2(z_p,z_a(i,:),@(XI,XJ) dist_fnct(XI,XJ,cov_inv(cs(i):ce(i),cs(i):ce(i))));
                %                 end
                
                cov_indxs = linspaceNDim(cs,ce,ce(1)+1-cs(1))';
                cov_indxs = cov_indxs(:);
                cov_part = cov_inv(cov_indxs,cov_indxs);
                logP(type_indx).mes = abs(pdist2Fast(z_p,z_a,cov_part));
            end
            
%             function D2 = dist_fnct(Xi,Xj,SIG)
%                 D2 = zeros(size(Xj,1),1);
%                 for row_indx=1:size(Xj,1)
%                    D2(row_indx) =  (Xi-Xj(row_indx,:))*SIG*(Xi-Xj(row_indx,:))';
%                 end  
%             end
            
        end
        
        function dz = makeDZ(obj,z,z_predict)
            %Compute measuremenr residuals for each pair of measurements of
            %the same type.
            
            %Get unique type of measurement input
            types = unique([z.type]);

            dz = SensorMeasurement.empty();
            for i=1:numel(types)
                mes_predicted = z_predict([z_predict.type] == types(i))
                mes_actual = z([z.type] == types(i));
                
                %We reuse the SensorMeasurment class here for convenience
                dz(i) = SensorMeasurement;
                dz(i).type = types(i);
                %The DZ matrix for a specific sensor type is 
                %as many rows as the measurement vector size
                %as many columns as there are actual measurements
                %as many 3rd dimension as there are predicted measurements
                dz(i).mes = zeros(numel(mes_actual(1).mes),numel(mes_actual),numel(mes_predicted));
                for j=1:numel(mes_predicted)
                   dz(i).mes(:,:,j) = [mes_actual.mes] - repmat(mes_predicted(j).mes,1,numel(mes_actual));
                end
            end
        end  
        
        
        function update_smt_matrix(obj, sensorId, newVal)
            findInd = find(obj.sensor_missing_thresh(:, 1) == sensorId);
            if ~isempty(findInd)
                obj.sensor_missing_thresh(findInd, 2) = newVal;
            end
        end
    end
    
    
    methods (Access=protected)
        
        function K = makeK(obj)
            %Make Optimal Kalman Gain 
            K = obj.P_predict*obj.H'/obj.inov_covariance;            %Optimal Kalman Gain
        end
        
        function P = makeP(obj)
            %Update covariance using Kalman Gain, we use the Joseph Form to
            %maintain PSD
            
            %Regular Form
            %P = obj.P_predict-obj.K*obj.H*obj.P_predict;
            
            %Joseph Form
            obs_noise = obj.MakeObservationNoise();
            obs_noise = obs_noise(obj.selectedcols,obj.selectedcols);
            I_m_KH = eye(size(obj.K,1)) - obj.K*obj.H;
            P = I_m_KH*obj.P_predict*I_m_KH' + obj.K*obs_noise*obj.K';
        end
    end
    
end
