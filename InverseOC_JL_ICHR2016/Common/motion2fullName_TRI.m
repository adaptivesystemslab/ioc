function motionSpecs = motion2fullName(motionAbbreviation)
    switch motionAbbreviation
        case 'ARTT_STD'
            motionFullName = 'tippy toes';
            kinematicChainDirection = 'AnkleToHip';
            
        case 'BKBR_SUP'
            motionFullName = 'back bridge';
            kinematicChainDirection = 'AnkleToHip';
            
        case 'GAIT_STD'
            motionFullName = 'gait';
            kinematicChainDirection = 'HipToAnkle_Mobile';
            
        case 'HAAO_STD'
            motionFullName = 'Hip abduction adduction, standing';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'HAAO_SUP'
            motionFullName = 'Hip abduction adduction, supine';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'HEFO_STD'
            motionFullName = 'Hip extension (leg moving backwards)';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'HFEO_STD'
            motionFullName = 'standing hip flexion (leg moving upwards)';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'HFEO_SUP'
            motionFullName = 'supine hip flexion (leg moving upwards)';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'KEFO_SIT'
            motionFullName = 'sitting leg extensions';
            kinematicChainDirection = 'HipToAnkle';

        case 'KEFO_SUP'
            motionFullName = 'supine leg extensions';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'KFEO_SIT'
            motionFullName = 'sitting leg flexion (backwards)';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'KFEO_STD'
            motionFullName = 'sitting leg flexion (backwards)';
            kinematicChainDirection = 'HipToAnkle';            
            
        case 'KFEO_SUP'
            motionFullName = 'supine knee downward motion';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'KHEF_STD'
            motionFullName = 'Knee-hip flexion (marching motion)';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'KHEF_SUP'
            motionFullName = 'supine knee/hip flexion';
            kinematicChainDirection = 'HipToAnkle';
            
        case 'LUNG_STD'
            motionFullName = 'standing lunges';
            kinematicChainDirection = 'AnkleToHip';
            
        case 'SQUA_STD'
            motionFullName = 'standing squats';
            kinematicChainDirection = 'AnkleToHip';
            
        case 'STAI_STD'
            motionFullName = 'stairs';
            kinematicChainDirection = 'HipToAnkle_Mobile';
  
        otherwise
            motionFullName = 'unknown motion';
            kinematicChainDirection = 'HipToAnkle';
            
            % unassigned: 
            dataType.list{11} = 'sitting sweep kick';
            dataType.list{12} = 'standing side kick';
            dataType.list{13} = 'standing circle kick';
            dataType.list{14} = 'standing front/side kick';
            dataType.list{15} = 'lying down diagonal leg raise';
            dataType.list{16} = 'lying down tracing out a circle';
            dataType.list{17} = 'lying down front/side kick';
            
            dataType.list{21} = 'walking';
            dataType.list{22} = 'stairs';
    end
    
    motionSpecs.motionFullName = motionFullName;
    motionSpecs.kinematicChainDirection = kinematicChainDirection;
    motionSpecs.initialPosition = motionAbbreviation(end-2:end);
end