function com_model = com_calc(model) %,q,dq,ddq)

% bodyFrameList = {model.bodies.name};
bodyMass = zeros(1,numel(model.bodies));
bodyCoM = zeros(3,numel(model.bodies));

for bodyNum = 1:numel(model.bodies)
    bodyMass(bodyNum) = model.bodies(bodyNum).m;
    bodyName = model.bodies(bodyNum).name;
    
    com = model.bodies(bodyNum).com;
    t = model.bodies(bodyNum).t;
    
    bodyCoM(:,bodyNum) = t(1:3,4) + t(1:3,1:3)*com;
end

com_model = (bodyCoM*bodyMass')./sum(bodyMass);

