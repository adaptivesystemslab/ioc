function [sigDof, stdDevFwd, stdDevCum] = sigDofSelectCumStdDev(templateVelocity, varDofkmeans)
    % Eliminate the low variance stuff by selecting 'sigDofs', based on
    % velocity ranges. The largest velocity range is selected as the sigDof
        
    varThesholdDof = 0.60;
    
    dof = size(templateVelocity, 2);
    stdDevFwd = zeros(dof, 1); %TODO should we switch over to velocity?
    
    for i = 1:dof
        signalInsp = templateVelocity(:, i);
        [muhat,sigmahat] = normfit(signalInsp);
        stdDevFwd(i) = sigmahat;
    end
    
    % pull the top 'thresholdVarDof' percentage for sigDof (currently not used)
    [stdDevSorted, stdDevMap] = sort(stdDevFwd, 'descend'); 
    stdDevSum = sum(stdDevSorted);
    stdDevPercent = stdDevSorted/stdDevSum;
    stdDevCum = cumsum(stdDevPercent);
    
    if varDofkmeans == 0
         sigDof = stdDevMap(1:find(stdDevCum > varThesholdDof, 1, 'first'));
    else
            % or how about if we used k-means clustering?
        [idx, c] = kmeans(stdDevFwd, varDofkmeans);
        [clusterMax, clusterInd] = max(c);
        sigDof = find(idx == clusterInd);
    end


%     stdDevName = jointNames(sigDof);
    
%     % if there are several sigDofs, use them all
%     sigDof = find(stdDevFwd > varThesholdDof);    
%     if isempty(sigDof)
%         % otherwise, they might all be below threshold...use just largest
%         [sigVal, sigDof] = max(stdDevFwd);
%     end
end