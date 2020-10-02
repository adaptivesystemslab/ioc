function testModelAllLinks(model, vis, linkTest)
    if ~exist('vis', 'var')
        vis = rlVisualizer('vis',640,480);
        vis.addModel(model);
    end

    if ~exist('linkTest', 'var')
        linkTest = 1:length(model.position);
    end
    
    ind = 0;
    for i = linkTest
        q = [(0):(pi/10):(pi) (pi):(-pi/10):(0)];
        for j = 1:length(q)
            model.position = zeros(size(model.position));
            model.position(i) = q(j);

%             applyMarkersToVisualization(vis, [], []);
            model.forwardPosition();
            model.forwardVelocity();
            model.forwardAcceleration();
            model.inverseDynamics();
            
            ind = ind + 1;
            q_array(ind, :) = model.position';
            tau_array(ind, :) = model.torque';
            
            vis.update();
            pause(0.1);
        end
    end
end