function sigDof = sigDofStdDev(velo, kmeansCount)
    dof = size(velo, 2);
    stdDevFwd = zeros(dof, 1); %TODO should we switch over to velocity?

    for i = 1:dof
        signalInsp = velo(:, i);
        [muhat,sigmahat] = normfit(signalInsp);
        stdDevFwd(i) = sigmahat;
    end

    % pull the top 'thresholdVarDof' percentage for sigDof (currently not used)
%     [stdDevSorted, stdDevMap] = sort(stdDevFwd, 'descend');
%     stdDevSum = sum(stdDevSorted);
%     stdDevPercent = stdDevSorted/stdDevSum;
%     stdDevCum = cumsum(stdDevPercent);
    %     sigDof = stdDevMap(1:find(stdDevCum > varThesholdDof, 1, 'first'));
    %     stdDevName = motionName(stdDevMap);

    % or how about if we used k-means clustering?
    [idx, c] = kmeans(stdDevFwd, kmeansCount);
    [clusterMax, clusterInd] = max(c);
    sigDof = find(idx == clusterInd);
end