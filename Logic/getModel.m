function model = getModel(trialInfo)
    switch trialInfo.baseModel
        case 'IIT'
            switch trialInfo.model
                case 'IIT_3DOF'
                    xmlPath = '../Libraries/rl/ik_framework/instance_iit/model/iit_v10_fixedbase_rightleg_sag.xml'; % sag 2D, 3dof half body model

                case 'IIT_7DOF'
                    xmlPath = '../Libraries/rl/ik_framework/instance_iit/model/iit_v10_fixedbase_right_sag.xml'; % sag 2D, 7dof half body model

                case 'IIT_17DOF'
                    xmlPath = '../Libraries/rl/ik_framework/instance_iit/model/iit_v10_fixedbase_right.xml'; % 3D, 17dof half body model
            end
            
            model = ModelRL_IIT(); 
            model.loadModel(xmlPath, trialInfo.path);
            
        case 'Jumping'
             switch trialInfo.model
                 case 'Jumping'
                     % load original IIT-based RL form
                     xmlPath = '../Libraries/rl/ik_framework/instance_jumping/model/JumpModel_IIT.xml';
                     model = ModelRL_Jumping_IIT();
                     model.loadModel(xmlPath, trialInfo.path);
                    
                 case 'Jumping2D'
                     % load Kevin's modified form
                    xmlPath = '../Libraries/rl/ik_framework/instance_jumping/model/JumpModel_Eul_inertia_2D.xml';
                    model = ModelRL_Jumping_KW(); 
                    model.loadModel(xmlPath, trialInfo);
             end

        case 'Expressive'
            switch trialInfo.model
                case 'ArmOnly'
                    xmlPath = '../Libraries/rl/ik_framework/instance_expressive/model/ioc_v4_rightarm_fixedbase.xml';
                    model = ModelRL_Expressive(trialInfo.model);
                    model.loadModel(xmlPath, trailInfo.path);
                    
                case 'BothArms'
                    xmlPath = '../Libraries/rl/ik_framework/instance_expressive/model/ioc_v4_upperbody.xml';
                    model = ModelRL_Expressive(trialInfo.model);
                    model.loadModel(xmlPath, trailInfo.path);
        case 'Healthy1'
            switch trialInfo.model
                case '4DOF'
                    xmlPath = '../Libraries/rl/ik_framework/instance_healthy1/model/healthy1_v3_rev1_sag.xml';
                    model = ModelRL_Healthy1();
                    model.loadModel(xmlPath, trialInfo.path);
            end
    end
end






