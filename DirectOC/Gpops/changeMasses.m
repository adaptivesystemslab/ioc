function outModel = changeMasses(inModel, masses)
    
    j = 1;
    for i=1:length(inModel.joints)
        currJoint = inModel.joints{i};
        % Check if link attached to i-th joint is not virtual
        if currJoint.mass ~= 0 && ~all(currJoint.com==0)
            currJoint.mass = masses(j);
            j = j+1;
        end
    end
    outModel = inModel;
end

