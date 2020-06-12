function [markerRmse, markerUsed, markerTotal, markerMav] = calc_rmse_mocap(featureSet_mocap, measurement_mask)
    rmseFct = setRmseFct();
    mavFct = setMavFct();
    
    markerRmse = [];
    markerMav = [];
    markerTotal = [];
    markerUsed = [];
    
    for ind_poserr = 1:length(featureSet_mocap.measurement_labels)
        %         ind = (1:3)+3*(2*ind_poserr-2);
        %         ind = (ind_poserr-1)*3+1:(ind_poserr)*3;
        ind = ind_poserr;
        currMask = measurement_mask(:, ind);
        currMesOrig = featureSet_mocap.measurement_input(:, ind).getMesArray;
        currMesEkf = featureSet_mocap.measurement_output(:, ind).getMesArray;

        %         % calculate euclidean norm
        %         currDiff = currMesEkf(:, 1:3) - currMesOrig(:, 1:3); % extract only position
        %         currDist = normVector(currDiff);
        %         currMasked = currDist(currMask(:, 1) == 1); % all 3 col of the mask should be the same
        %
        %         % calculate RMSE
        %         currRmse = sqrt(sum(currMasked.^2)/(numel(currMasked)));

        currRmse = rmseFct(currMesEkf(currMask == 1, 1:3), currMesOrig(currMask == 1, 1:3));
        currMav = mavFct(currMesEkf(currMask == 1, 1:3), currMesOrig(currMask == 1, 1:3));
        
        if isnan(currRmse)
            currRmse = 0;
            currMav = 0;
        end
        
        markerRmse(ind_poserr) = currRmse;
        markerMav(ind_poserr) = currMav;  
        markerTotal(ind_poserr) = size(currMesOrig, 1);
        markerUsed(ind_poserr) = sum(currMask == 1);
    end
end