function feature_full = feature_collate(traj_load, feature_load_names)
    % expects an input struct that contains a preset of features. if the
    % feature is not passed in, then assign a blank array as the value
    p = inputParser;
    p.KeepUnmatched = true;
    
    % set up a blank array for data that has not been loaded
    blankArray = zeros(size(traj_load.q)); % assuming 'q' always will exist
    
    for ii = 1:length(feature_load_names)
        addOptional(p, feature_load_names{ii}, blankArray);
    end

    parse(p, traj_load); % perform the file checking
    feature_full = p.Results;    
end