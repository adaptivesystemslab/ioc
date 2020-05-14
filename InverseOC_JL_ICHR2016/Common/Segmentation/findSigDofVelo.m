function [sigVelo, sigDofInd] = findSigDofVelo(velo)
    % find the most sig dof based on max var
    
%     velo = jointAngleData.jointVelo';
    veloVar = var(velo);
    [maxInd, sigDofInd] = max(veloVar);
%     sigDof = sigOfConsideration(sigDofInd);
    sigVelo = [velo(:, sigDofInd)]';
end