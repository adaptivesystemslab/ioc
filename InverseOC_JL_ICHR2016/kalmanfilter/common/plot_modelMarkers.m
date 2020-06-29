function plotModelMarkers(obj)
    vis = rlVisualizer('vis',640,480);
    vis.addModel(obj.model);

    applyMarkersToVisualization(vis, obj.measurement_label, obj.measurement_indArray, 1);

    obj.model.forwardPosition();
    vis.update();
    pause;
end