function sigDof = segmentModiftZvcMethod_preamble(velo)

sigDof = 0;

if ~isempty(velo)

    % peak ZVC method
    %         windowMode = 1;
    veloMode = 'sigdof'; % norm sigdof
    
    switch veloMode
        case 'norm'
            % calculate the norm velocity. use all the dofs
            sigDof = 1:size(velo, 2);
            
        case 'sigdof'
            for indVelo = 2:size(velo, 2)
                % try different kmeans count to determine a good number
                % for sigdof count
                try
                    [sigDof, stdDevFwd, stdDevCum] = sigDofSelectCumStdDev(velo, indVelo);
                    
                    if length(sigDof) < 1
                        % don't keep more than 4 sigdofs
                        break
                    end
                catch err
                    
                end
            end
    end
end
