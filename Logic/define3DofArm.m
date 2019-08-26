function arm = define3DofArm(modelType, masses)

    if ~exist('modelType', 'var')
        modelType = 'rt';
    end
    
    if ~exist('masses', 'var')
        masses = [0, 1, 1];
    end
    
    addpath(genpath('../Utils'));
       
    tensor = struct;
    tensor.type = 'cylinder';
    tensor.direction = 'x';
    tensor.length = 1;
    tensor.radius = 0.4;
    tensor.translate = true;
    
    chain = struct;
    chain(1).a = 0;
    chain(1).d = 0;
    chain(1).r = 0;
    chain(1).theta = 0;
    chain(1).virtual = 1;
    chain(1).endEffector = 0;
    
    chain(2).a = pi/2;
    chain(2).d = 0;
    chain(2).r = 0;
    chain(2).theta = 0;
    chain(2).virtual = 0;
    chain(2).m = masses(2);
    chain(2).com = [.5 0 0];
    chain(2).tensor = tensor;
    chain(2).endEffector = 0;
    
    chain(3).a = 0;
    chain(3).d = 0;
    chain(3).r = 1;
    chain(3).theta = 0;
    chain(3).virtual = 0;
    chain(3).m = masses(3);
    chain(3).com = [.5 0 0];
    chain(3).tensor = tensor;
    chain(3).endEffector = 0;
                  
    gravity = [0 0 -9.81];
    
   switch lower(modelType)
       case 'symoro'
           armModel = ArmModelSymoro(chain, gravity);
       case 'rt'
           armModel = ArmModelRT(chain);
       case 'rl'
           armModel = ArmModelRL();
           armModel.createModel(chain);
           armModel.initModel();
    end 
    
    arm = armModel;
end