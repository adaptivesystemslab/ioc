function model = getModel(trialInfo)
%     if ~exist('masses', 'var')
%         masses = [0, 1, 1];
%     end
%     
    switch trialInfo.baseModel
        case 'IIT'
            model = ArmModelRL();
            model.loadAndSetupIIT(trialInfo.path, trialInfo);
            model.initModel();
            
        case 'Jumping'
            model = ArmModelRL();
            model.loadAndSetupJumping(trialInfo.path, trialInfo);
            model.initModel();
    end
end

