function testJointRange(modelInstance)
    mdl = modelInstance.model;
    
    vis = rlVisualizer('vis',640,480);
    mdl.forwardPosition();
    vis.addModel(mdl);
    vis.update();
    
%     jointRangeOfMotion = modelInstance.jointRangeOfMotion;

    
    
    
    
    
     jointRangeOfMotion = {...
           {'joint_rankle_2',      [  -18, 33]} ...
         {'joint_lankle_2',       [ -33, 18]} ...
        };

    skips = 3;
    render(vis, mdl, jointRangeOfMotion, skips);
    
    
    
end

function render(vis, mdl, jointRangeOfMotion, skips)
    allJointStr = {mdl.joints.name};
    
    for i = 1:length(jointRangeOfMotion)
        currJointLimArray = jointRangeOfMotion{i};
        currQName = currJointLimArray{1};
        currQRan = currJointLimArray{2};
        
%         q = deg2rad([(0):(skips):(currQRan(1)) (currQRan(1)):(-skips):(0) ...
%                      (0):(skips):(currQRan(2)) (currQRan(2)):(-skips):(0)]);
        q = deg2rad([(currQRan(1)):(skips):(currQRan(2)) (currQRan(2)):(-skips):(currQRan(1))]);
        
        indx = find(ismember(allJointStr,currQName),1);
        
        for j = 1:length(q)
            mdl.position(4:end) = 0;
            mdl.position(indx) = q(j);
            
            mdl.forwardPosition();
            vis.update();
            pause(0.1);
        end
    end
end
