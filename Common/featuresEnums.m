classdef featuresEnums
    % to add a new cost function, add it here in the enum, then add a
    % corresponding entry in featureCalc.featureCalc();
    
    % also implement: if no weights are given, load the masses of the frames/joints
    % if evenlyl distributed mass if not COM-based things
    % also, list of all frames and joints to default into
    
    enumeration
        cartVeloSumSqu % sumsq cart velo of all frames (or some subset)
        cartAccelSumSqu % sumsq cart accel of all frames (or some subset)
        cartJerkSumSqu % sumsq cart jerk of all frames (or some subset)
        cartCurvatureSumSqu % sumsq curvature of all frames (or some subset)
        cartRadCurvatureSumSqu % inverse of cartCurvatureSumSqu 
        cartBoundingBoxSumSqu 
        cartBoundingVolumeSumSqu 
        cartDisplacementSumSqu 
        cartQuantityMotionSumSqu % weighted end eff velo
        cartWeightEffortSumSqu 
        cartTimeEffortSumSqu 
        cartSpaceEffortSumSqu 
        cartFlowEffortSumSqu 
        
        centreMassSumSqu % overall centre of mass
        centreMassVeloSumSqu % overall centre of mass velo
        centreMassDisplacementSumSqu 
        
        angVeloSumSqu
        angAccelSumSqu
        angJerkSumSqu
        angCurvatureSumSqu
        angRadCurvatureSumSqu
        angQuantityMotionSumSqu 
        angWeightEffortSumSqu 
        
        torqueSumSqu
        torqueVeloSumSqu
        torqueAccelSumSqu
        
        kineticEnergySumSqu
        angPowerSumSqu
        
        extensivenessMaxSumSqu 
        extensivenessSumSumSqu 
        
        shapeDirectionSumSqu 
        
        cartDistToTarget
        rotDistToTarget
        
        
        
       
        cartVeloAxis % cart velo of a given frame, in either xyz
        centreMassVeloRelativeToFrame % com velo relative to some cart frame velo
    end
end