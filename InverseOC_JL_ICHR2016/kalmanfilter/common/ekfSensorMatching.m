function [ekf_markerMask, ekf_eventhandler, ekf_markerTally] = ekfSensorMatching(featureSet)
    ekf_match = featureSet.ekf_match;
    matchMatrix = featureSet.measurement_ekf_match;
    output_matching = featureSet.measurement_output_match;

    %             ekf_labels = obj.measurement_labels;
    %             all_labels = dataInstance.measurement_labels;

    % calculate measurement mask if markers are swapped/missing
    % ekf_match: len is mdl.sensors, vals is indices of meas.allSens
    % matchMatrix: len is mdl.sensors, vals is indices of meas.allSen
    % outputMatch: len is mdl.allSens, vals is indices of mdl.sensor
    lenTime = length(featureSet.time);
    arrayTime = (1:lenTime)';
    maskNone = ones(lenTime, length(output_matching));
    maskCorrect = zeros(lenTime, length(output_matching));
    maskSwapMissing = zeros(lenTime, length(output_matching));
    maskSwapped = zeros(lenTime, length(output_matching));
    maskMissing = zeros(lenTime, length(output_matching));
    
    markersTotal = zeros(1, length(output_matching)); % tally number of markers/frame total
    markersCorrect = zeros(1, length(output_matching)); % tally number of markers/frame that are correct
    markersSwapped = zeros(1, length(output_matching)); % tally number of markers/frame that are missing (ie = 0)
    markersMissing = zeros(1, length(output_matching)); % tally number of markers/frame that are swapped (ie = different number than matchmatrix)
    
    for i = 1:length(output_matching)
        currSenInd = output_matching(i);
        if currSenInd ~= 0
            currMes = matchMatrix(output_matching(i));
            currEkfMatch = ekf_match(:, currSenInd);
            
            markersTotal(i) = length(currEkfMatch);
             
            entries_c = find(currEkfMatch == currMes);
            markersCorrect(i) = length(entries_c);
            maskCorrect(entries_c, i) = 1;
            
            entries_sm = find(currEkfMatch ~= currMes);
            maskSwapMissing(setxor(arrayTime, entries_sm), i) = 1; % maintained in the same order as ekf sensor (orig)
            
            entries_m = find(currEkfMatch(entries_sm) == 0);
            markersMissing(i) = length(entries_m);
            maskMissing(setxor(arrayTime, entries_sm(entries_m)), i) = 1;

%             inds = find(currEkfMatch > 0);
            entries_s = find(currEkfMatch(entries_sm) > 0);
            markersSwapped(i) = length(entries_s);
            maskSwapped(setxor(arrayTime, entries_sm(entries_s)), i) = 1;
        end
        
        
    end
    
    ekf_eventhandler = featureSet.ekf_eventhandler;
    ekf_markerMask.none = maskNone;
    ekf_markerMask.swappedmissing = maskSwapMissing;
    ekf_markerMask.missing = maskMissing;
    ekf_markerMask.swapped = maskSwapped;
    ekf_markerMask.correct = maskCorrect;
    ekf_markerTally.total = markersTotal;
    ekf_markerTally.correct = markersCorrect;
    ekf_markerTally.swapped = markersSwapped;
    ekf_markerTally.missing = markersMissing;
end

    %         if currMes > size(ekf_match, 2) || currMes == 0
%             continue
%         end

    

    % tally number of ekf events
