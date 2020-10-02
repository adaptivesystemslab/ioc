function [observation, fhmmIntermediate] = updateDataStruct(observation, observationNew, sysparam, fhmmIntermediate)   

    obsTime = observation.obsTime;
    obsData = observation.obsData;
    obsVelo = observation.obsVelo;
    
    obsTimeNew = observationNew.obsTime';
    obsDataNew = observationNew.obsData';
    obsVeloNew = observationNew.obsVelo';

    if observation.currDataLength == 0
        % first entry, so need some padding to wind up the filter
        sysparam.filter.filterObservation.padLength = 50;
        
%         [newData, zi_ang] =  filter_butterworth_function(obsDataNew', sysparam.filter.filterType, ...
%             sysparam.filter.typeObservation, sysparam.filter.bwFreq, sysparam.filter.bwOrder, [], 50);
%         [newVelo, zi_velo] =  filter_butterworth_function(obsVeloNew', sysparam.filter.filterType, ...
%             sysparam.filter.typeObservation, sysparam.filter.bwFreq, sysparam.filter.bwOrder, [], 50);

        [newData] = obsDataNew';
        [newVelo] = obsVeloNew';
        zi_ang = 0;
        zi_velo = 0;
        
        
%         [newData, zi_ang] = sysparam.filter.filterObservation.apply(observationNewData');
%         [newVelo, zi_velo] = sysparam.filter.filterObservation.apply(observationNewDataExtra');
        fhmmIntermediate.zi_ang = zi_ang;
        fhmmIntermediate.zi_velo = zi_velo;
        
        % initalizing the array with new data
        obsTime = obsTimeNew;
        obsData = newData';
        obsVelo = newVelo';
        
        observation.currDataLength = length(obsTime);
        
        veloInd = 1:length(obsTimeNew);
        newVeloInd = 1:size(obsData, 2);
        indDataAdding = newVeloInd;
        %     veloInd = 1:observation.currDataLength;
    else
%         sysparam.filter.filterObservation.padLength = 0;

% % %         [newData, zi_ang] =  filter_butterworth_function(observationNewData', sysparam.filter.filterType, ...
% % %             sysparam.filter.typeObservation, sysparam.filter.bwFreq, sysparam.filter.bwOrder, fhmmIntermediate.zi_ang, 0);
% % %         [newVelo, zi_velo] =  filter_butterworth_function(observationNewDataExtra', sysparam.filter.filterType, ...
% % %             sysparam.filter.typeObservation, sysparam.filter.bwFreq, sysparam.filter.bwOrder, fhmmIntermediate.zi_velo, 0);
%         [newData, zi_ang] = sysparam.filter.filterObservation.apply(observationNewData', fhmmIntermediate.zi_ang);
%         [newVelo, zi_velo] = sysparam.filter.filterObservation.apply(observationNewDataExtra', fhmmIntermediate.zi_velo);

        [newData] = obsDataNew';
        [newVelo] = obsVeloNew';
        zi_ang = 0;
        zi_velo = 0;

        fhmmIntermediate.zi_ang = zi_ang;
        fhmmIntermediate.zi_velo = zi_velo;
        
        % add the new data into the array
        lenDataAdding = length(obsTimeNew);
        indDataAdding = observation.currDataLength + 1 : observation.currDataLength + lenDataAdding;
        
        % in case the matrix is not big enough...
        obsTime(indDataAdding(end)) = 0;
        obsData(1, indDataAdding(end)) = 0;
        obsVelo(1, indDataAdding(end)) = 0;
         
        obsTime(indDataAdding) = obsTimeNew;
        obsData(:, indDataAdding) = newData';
        obsVelo(:, indDataAdding) = newVelo';
       
        observation.currDataLength = observation.currDataLength + lenDataAdding;
        
        veloInd = [indDataAdding(1)-1 indDataAdding];
        newVeloInd = 2:length(veloInd);
%     veloInd = 1:observation.currDataLength;
    end
    
%     if observation.currDataLength >= length(obsTime)
%         % double the length of the data array since it is now full
%         indToExtendTo = observation.currDataLength * 2;
%         
%         obsTime(indToExtendTo) = 0;
%         obsData(1, indToExtendTo) = 0;
%         obsVelo(1, indToExtendTo) = 0;
%     end
    
%     if strcmpi(sysparam.filter.typeObservation, 'movingAverage')
%         newVelo = veloCalc(obsTime(veloInd)', obsData(observation.posIndex, veloInd)', sysparam.filter.typeObservation, sysparam.filter)';
%         obsVelo(:, indDataAdding) = newVelo(:, newVeloInd);
%     else
%         veloInd = 1:observation.currDataLength; % override
%         obsVelo = veloCalc(obsTime(veloInd)', obsData(observation.posIndex, veloInd)', sysparam.filter.typeObservation, sysparam.filter)';
%     end
    
    observation.obsTime = obsTime;
    observation.obsData = obsData;
    observation.obsVelo = obsVelo;
end
