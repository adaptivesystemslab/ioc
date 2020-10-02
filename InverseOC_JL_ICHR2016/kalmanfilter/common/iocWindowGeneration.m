function indToUse_iocWindow = iocWindowGeneration(iocRec, featureSetIocFull, segmentInfoDoc)
    switch iocRec.windowMethod
        case 'sliding'
            ind = 0;
            for ind_windowing = 1:iocRec.windowAdvance:numel(featureSetIocFull.time)
                ind = ind + 1;
                startInd = ind_windowing;
                endInd = ind_windowing + iocRec.windowLength - 1;

                if endInd <= numel(featureSetIocFull.time)
                    indToUse_iocWindow(ind, 1) = startInd;
                    indToUse_iocWindow(ind, 2) = endInd;
                end
            end

        case 'segments'
            for ind_windowing = 1:numel(segmentInfoDoc)
                [startVal, startInd] = findClosestValue(segmentInfoDoc(ind_windowing).timeStart, featureSetIocFull.time);
                [endVal, endInd] = findClosestValue(segmentInfoDoc(ind_windowing).timeEnd, featureSetIocFull.time);

                indToUse_iocWindow(ind_windowing, 1) = startInd;
                indToUse_iocWindow(ind_windowing, 2) = endInd;
            end
    end