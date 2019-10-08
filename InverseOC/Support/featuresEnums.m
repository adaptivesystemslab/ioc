classdef featuresEnums
    % to add a new cost function, add it here in the enum, then add a
    % corresponding entry in featureCalc.featureCalc();
    
    % also implement: if no weights are given, load the masses of the frames/joints
    % if evenlyl distributed mass if not COM-based things
    % also, list of all frames and joints to default into
    
    enumeration
        cartVeloSumSqu 
        cartAccelSumSqu
        cartJerkSumSqu
        cartCurvatureSumSqu
        cartRadCurvatureSumSqu
        cartBoundingBoxSumSqu
        cartBoundingVolumeSumSqu 
        cartDisplacementSumSqu 
        cartQuantityMotionSumSqu 
        cartWeightEffortSumSqu 
        cartTimeEffortSumSqu 
        cartSpaceEffortSumSqu 
        cartFlowEffortSumSqu 
        
        centreMassSumSqu 
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
    end
end