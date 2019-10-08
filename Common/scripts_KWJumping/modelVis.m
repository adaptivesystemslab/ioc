mdl_eul = rlCModel('JumpModel_Eul.xml');

vis = rlVisualizer('vis', 480, 960);
vis.addModel(mdl_eul);
mdl_eul.forwardPosition;
vis.update;

clear vis;