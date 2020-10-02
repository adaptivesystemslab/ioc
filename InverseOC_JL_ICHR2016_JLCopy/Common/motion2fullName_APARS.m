function motionSpecs = motion2fullName(motionAbbreviation)
    switch motionAbbreviation
        case 'KEFO_SIT'
            motionFullName = 'Sitting, leg extensions';
            kinematicChainDirection = 'HipToAnkle';
            description = 'Knee extension/flexion motion, extending from 90d to 180d.';
            
        case 'STSO_SIT'
            motionFullName = 'Sitting, sit-to-stand';
            kinematicChainDirection = 'HipToAnkle';
            description = '(As natural as possible)';
            
        case 'HARO_SIT'
            motionFullName = 'Sitting, axial rotation'; % low sample count
            kinematicChainDirection = 'HipToAnkle';
            description = 'Rotate leg inward as far as possible, then back.';
            
        case 'LUNG_STD'
            motionFullName = 'Lunges';
            kinematicChainDirection = 'AnkleToHip';
            description = 'Right leg step forward and lean forward until right knee is at 90d, while maintaining a straight back. Return to resting position by pushing off with right leg.';
            
        case 'SQUA_STD'
            motionFullName = 'Squats';
            kinematicChainDirection = 'AnkleToHip';
            description = 'Keeping upper body upright, bend knees so that the bent of the knee is in the same direction as the feet. Try to keep the lower leg as straight as possible and bend mainly in the hip.  ';
            
        case 'MARC_STD'
            motionFullName = 'Standing, knee raise';
            kinematicChainDirection = 'HipToAnkle';
            description = 'marching motion, with knee coming up to the waist.';
            
        case 'KICK_STD'
            motionFullName = 'Standing, front kicks'; % low sample count
            kinematicChainDirection = 'HipToAnkle';
            description = '';
            
        case 'KHEF_SUP'
            motionFullName = 'Lying down, knee/hip flexion';
            kinematicChainDirection = 'HipToAnkle';
            description = 'This is best done laying on your back essentially sliding your heel along the bed towards your buttocks.';
            
        case 'HFEO_SUP'
            motionFullName = 'Lying down, straight leg raises';
            kinematicChainDirection = 'HipToAnkle';
            description = 'Lay on your back with non-operative lower limb bent and foot resting on bed. This acts to give stability against spinal rotation. Then, keeping the operative knee as straight as possible, flex at the hip to raise the entire leg up off the bed about 1-1.5 inches. Goal: keep knee extended fully even as you lower the leg down.';
            
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
    
    motionSpecs.description = description;
    motionSpecs.motionFullName = motionFullName;
    motionSpecs.kinematicChainDirection = kinematicChainDirection;
    motionSpecs.initialPosition = motionAbbreviation(end-2:end);
end