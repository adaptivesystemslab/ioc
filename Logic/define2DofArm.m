function arm = define2DofArm(modelType)
    
    if ~exist('modelType', 'var')
        modelType = 'rt';
    end
    
    tensor = struct;
    tensor.type = 'default';
    tensor.translate = true;
    
    chain = struct;
    chain(1).a = pi/2;
    chain(1).d = 0;
    chain(1).r = 0;
    chain(1).theta = 0;
    chain(1).virtual = 0;
    chain(1).m = 1;
    chain(1).com = [.5 0 0];
    chain(1).tensor = tensor;
    chain(1).endEffector = 0;
    
    chain(2).a = 0;
    chain(2).d = 0;
    chain(2).r = 1;
    chain(2).theta = 0;
    chain(2).virtual = 0;
    chain(2).m = 1;
    chain(2).com = [.5 0 0];
    chain(2).tensor = tensor;
    chain(2).endEffector = 0;
    
    gravity = [0 0 -9.81];
    
    if strcmp(lower(modelType), 'symoro')
            armModel = ArmModelSymoro(chain, gravity);
    else
            armModel = ArmModelRT(chain);
    end
    arm = armModel;

end

