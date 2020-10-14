function com_visualizer_test(model,vis)

% vis = rlVisualizer('vis', 640, 960);
% vis.addModel(model);
% vis.update;


% figure(1); clf; hold on; grid on;

for bodyNum = 1:numel(model.bodies)
    com = model.bodies(bodyNum).com;
    t = model.bodies(bodyNum).t;
    com_global = t(1:3,4) + t(1:3,1:3)*com;
%     com_global = t(1:3,4) + com;
    vis.addMarker(model.bodies(bodyNum).name,com_global);
end
% 
% clear vis;
