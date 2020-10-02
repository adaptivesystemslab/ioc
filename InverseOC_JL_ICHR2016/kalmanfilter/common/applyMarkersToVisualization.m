function applyMarkersToVisualization(vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    % add markers
    colours.mdl                  = [0 0 0 0];     % model sensor name
    colours.mesEkfExactMatch     = [0.5 0 0.5 1];  % in mes and assigned to an ekf entry
    colours.mesEkfNotExactMatch  = [0 1 1 1];  % in mes and assigned to an ekf entry
    colours.mesNotEkf            = [1 1 0 0.5];     % in mes but not assigned to an ekf entry
    colours.trc                  = [1 0 0 0.1];     % in trcdata but not in mes
    colours.axis                 = [0 0 0 0.2];    % axis labels
    
    if ~exist('cur_ass_mes', 'var')
        cur_ass_mes = [];
    end
    
    if ~exist('matchMatrix', 'var')
        matchMatrix = [];
    end
    
    % add coordinate axes
    plotAxisMarker(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix);
    
    if ~isempty(mdl)
        % visualize model markers
        plotModelMarker(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix);
        plotModelJoints(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix);
    end
    
    if ~isempty(dataInstance)
        if ~isempty(dataInstance.data)
            if isempty(cur_ass_mes) && isempty(matchMatrix)
                colours.trc                  = [1 0 0 1];
                plotTrcMarkerAll(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix);
            else
                % visualize the unused markers
                plotTrcMarkerRestrict(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix);
                
                % visualize all markers in dataInstance.measurement: ekf vs not
                plotEkfMarker(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
            end
        else
            % if dataInstance.data is empty, check the measurement array for any
            % position data and plot that
            plotMeasurementPosition(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix);
        end
    end
end

function plotAxisMarker(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    vis.addMarker('x-axis', [1 0 0], colours.axis);
    vis.addMarker('y-axis', [0 1 0], colours.axis);
    vis.addMarker('z-axis', [0 0 1], colours.axis);
end

function plotModelMarker(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    for j = 1:length(mdl.sensors)
        name = ['mSensor_' mdl.sensors(j).name];
        pos = mdl.sensors(j).transform(1:3, 4);
        vis.addMarker(name, pos, colours.mdl);
    end
end

function plotModelJoints(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    for j = 1:length(mdl.joints)
        name = ['mJoint_' mdl.joints(j).name];
        pos = mdl.joints(j).frame_in.t(1:3, 4);
        vis.addMarker(name, pos, colours.mdl);
    end
end

function plotTrcMarkerAll(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    trcFields = fieldnames(dataInstance.data);
    for m = 1:length(trcFields)
        name = ['trc_' trcFields{m}];
        pos = dataInstance.data.(trcFields{m})(frameInd, :);
        vis.addMarker(name, pos, colours.trc);
    end
end

function plotMeasurementPosition(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    for m = 1:length(dataInstance.measurement_labels)
        name = ['trc_' dataInstance.measurement_labels{m}];
        if dataInstance.measurement(1, m).type == 1
            pos = dataInstance.measurement(frameInd, m).getMesArray;
            vis.addMarker(name, pos, colours.mesEkfExactMatch);
        end
    end
end

function plotTrcMarkerRestrict(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    if length(dataInstance.data) > 1
        % if there is an array of data, it is probably imu data
        return
    end

    trcFields = fieldnames(dataInstance.data);
    for m = 1:length(trcFields)
        % if it's in the main sensor array, don't plot it here
        sensorNameSplit = strsplit(trcFields{m}, '_');
        if sum(strcmpi(dataInstance.measurement_labels, trcFields{m}) == 1)
            % if it's already in measurement, don't plot it
            continue;
        elseif strcmpi(sensorNameSplit{end}, 'JC')
            % if it's synthetic, don't plot it
            continue;
        end
        
        name = ['trc_' trcFields{m}];
        pos = dataInstance.data.(trcFields{m})(frameInd, :);
        vis.addMarker(name, pos, colours.trc);
    end
end

function plotEkfMarker(colours, vis, mdl, dataInstance, frameInd, cur_ass_mes, matchMatrix)
    lenMesLabels = 1:length(dataInstance.measurement_labels);
    for j = lenMesLabels
        mesInd = matchMatrix(j, 2);
        ekf_ind = find(mesInd == cur_ass_mes);
        
        if dataInstance.measurement(frameInd, mesInd).type ~= 1
            continue
        end
            
        pos = dataInstance.measurement(frameInd, mesInd).mes(1:3);

        if ~isempty(ekf_ind)
            mdlSensName = mdl.sensors(ekf_ind).name;
        else
            mdlSensName = '';
        end
        mesSensName = dataInstance.measurement_labels{mesInd};
        name = ['mes_' mesSensName];

        if ~isempty(ekf_ind) && strcmpi(mdlSensName, mesSensName)
            % this measurement entry is matched to an ekf entry
            %             name = ['mes_' dataInstance.measurement_labels{j} '_' 'ekf_' mdl.sensors(ekf_ind).name];
            vis.addMarker(name, pos, colours.mesEkfExactMatch);
        elseif ~isempty(ekf_ind)
            vis.addMarker(name, pos, colours.mesEkfNotExactMatch);
        else
            vis.addMarker(name, pos, colours.mesNotEkf);
        end
    end
end