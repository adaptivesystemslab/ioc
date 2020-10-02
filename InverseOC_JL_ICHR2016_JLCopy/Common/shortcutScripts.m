function shortcutScripts(fctcall)
    % A collection of shortcuts used in the MATLAB shortcut bar. The more
    % complicated functions will be placed here, whereas if it's just one
    % line or something, I just left it in the shortcut bar
    
    switch fctcall
        case 'AlignFig'
            % In dual-monitor setups, sometimes the figures are displayed
            % offscreen. A recommended fix:            
            AlignDefaultFig
            
        case 'MLintClear'
            % Sometimes MLint would just fail for no good reason. This will
            % clear the MLintFailureFile to reset that            
            mlintClear
            
        otherwise
            fprintf('SCS: Function not found. Check call code\n');
    end
end