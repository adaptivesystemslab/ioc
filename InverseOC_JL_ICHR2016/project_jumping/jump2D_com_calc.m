function com_model = jump2D_com_calc(model) %,q,dq,ddq)

visu = 0;

% model = inputData(model, q, dq, ddq);
% model.forwardPosition();
% model.forwardVelocity();
% model.forwardAcceleration();


if(visu)
    vis = rlVisualizer('vis', 640, 960);
    vis.addModel(model);
    vis.update;
end

bodyFrameList = {model.bodies.name};
bodyMass = zeros(1,numel(model.bodies));
bodyCoM = zeros(3,numel(model.bodies));

for bodyNum = 1:numel(model.bodies)
    bodyMass(bodyNum) = model.bodies(bodyNum).m;
    bodyName = model.bodies(bodyNum).name;
    
    com = model.bodies(bodyNum).com;
    t = model.bodies(bodyNum).t;
    
    % Specific cases
    if(contains(bodyName,'back')) % back5 rotated up 90 degrees
        bodyCoM(:,bodyNum) = t(1:3,4) + roty(-pi/2)*t(1:3,1:3)*com;
    elseif(contains(bodyName,'rankle'))
        bodyCoM(:,bodyNum) = t(1:3,4) + rotx(pi)*roty(pi)*t(1:3,1:3)*rotz(pi/2)*roty(pi)*com; % not efficient but works
    elseif(contains(bodyName,'lankle'))
        bodyCoM(:,bodyNum) = t(1:3,4) + roty(pi)*t(1:3,1:3)*rotz(pi/2)*com;
    elseif(contains(bodyName,'relbow')) % right side bodies flipped in Y axis
        bodyCoM(:,bodyNum) = t(1:3,4) + roty(pi)*t(1:3,1:3)*rotz(-pi/2)*com;
    % rest of right side bodies flipped in Y axis
    elseif(bodyName(1) == 'r')
        bodyCoM(:,bodyNum) = t(1:3,4) + roty(pi)*t(1:3,1:3)*rotz(pi)*com;
    else
        bodyCoM(:,bodyNum) = t(1:3,4) + t(1:3,1:3)*com;
    end
    
    if(visu)
        vis.addMarker(model.bodies(bodyNum).name,bodyCoM(:,bodyNum));
    end
end

com_model = (bodyCoM*bodyMass')./sum(bodyMass);

if(visu)
    vis.addMarker('model_CoM',com_model);
%     clear vis;
end
