function arm = define4DofArm(modelType)

    if ~exist('modelType', 'var')
        modelType = 'symoro';
    end
    
    tensor = struct;
    tensor.type = 'cylinder';
    tensor.direction = 'z';
    tensor.length = 1;
    tensor.radius = 0.4;
    tensor.translate = true;
    
    chain = struct;
    chain(1).a = pi/2;
    chain(1).r = 0;
    chain(1).theta = pi/2;
    chain(1).d = 0;
    chain(1).virtual = 1;
    chain(1).endEffector = 0;
    
    chain(2).a = pi/2;
    chain(2).r = 0;
    chain(2).theta = -pi/2;
    chain(2).d = 0;    
    chain(2).virtual = 1;
    chain(2).endEffector = 0;

    chain(3).a = -pi/2;
    chain(3).r = 0;
    chain(3).theta = -pi/2;
    chain(3).d = -1;
    chain(3).m = 2;
    chain(3).com = [0 0 0.5];
    chain(3).virtual = 0;
    chain(3).endEffector = 0;
    
    chain(4).a = pi/2;
    chain(4).r = 0;
    chain(4).theta = 0;
    chain(4).d = 0;
    chain(4).virtual = 1;
    chain(4).endEffector = 0;
    
    chain(5).a = -pi/2;
    chain(5).r = 0;
    chain(5).theta = 0;
    chain(5).d = -1;             
    chain(5).virtual = 0;
    chain(5).m = 1.5;
    chain(5).com = [0 0 0.5];
    chain(5).endEffector = 1;
    
    gravity = [0 0 -9.81];
    
    switch lower(modelType)
        case 'symoro'
            chain(3).tensor = tensor;
            chain(5).tensor = tensor;
            armModel = ArmModelSymoro(chain, gravity);
            
        case 'rt'
            tensor.translate = false;
            chain(3).tensor = tensor;
            chain(5).tensor = tensor;
            armModel = ArmModelRT(chain);
            
        case 'rl'
            tensor.translate = false;
            chain(3).tensor = tensor;
            chain(5).tensor = tensor;
            armModel = ArmModelRL();
            armModel.createModel(chain);
            armModel.initModel();
    end        

    arm = armModel;
end