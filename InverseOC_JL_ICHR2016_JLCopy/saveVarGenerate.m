% compress all the variables needed for plotting into savevar

saveVar.run_mode = run_mode;
saveVar.indToUse_window = indToUse_window;
saveVar.feature_full = feature_full;
saveVar.ccost_full = ccost_full;
saveVar.output_inverse = output_inverse;
saveVar.t_recon_plot_array = t_recon_plot_array;
saveVar.elapsedTime = elapsedTime;
saveVar.q_recon_plot_array = q_recon_plot_array;
saveVar.dq_recon_plot_array = dq_recon_plot_array;
saveVar.ddq_recon_plot_array = ddq_recon_plot_array;
saveVar.minRmseInd_plot = minRmseInd_plot;
saveVar.feature_recon_local = feature_recon_local;
saveVar.feature_win_save = feature_win_save;
saveVar.minRmseIndArray = minRmseIndArray;
saveVar.segmentInfo = segmentInfo;
saveVar.cost_function_names = cost_function_names;
saveVar.runSettings = runSettings;
saveVar.lenToUseT = lenToUseT;
saveVar.currFilestack = currFilestack;
 

paramSave.dataset = param.dataset;
paramSave.segment_only_windows = param.segment_only_windows;
paramSave.win_shift = param.win_shift;
paramSave.win_length = param.win_length;
paramSave.const_x = param.const_x;
paramSave.const_y = param.const_y;
paramSave.t_full = param.t_full;
paramSave.t_spline = param.t_spline;
paramSave.x_knot_set = param.x_knot_set;
% paramSave.spline_fit = param.spline_fit;

paramSave.spline_length = param.spline_length;
paramSave.half_phase_construction = param.half_phase_construction;
paramSave.doc_sim_sp_dq_endcond = param.doc_sim_sp_dq_endcond;
paramSave.doc_sim_cf_dq_constraint = param.doc_sim_cf_dq_constraint;
paramSave.doc_sim_win_length = param.doc_sim_win_length;
paramSave.ioc = param.ioc;
paramSave.doc_pivot = param.doc_pivot;
paramSave.corrThreshold = param.corrThreshold;
paramSave.ioc_add_knot_as_cf_constraints = param.ioc_add_knot_as_cf_constraints;
paramSave.ioc_sp_dq_endcond = param.ioc_sp_dq_endcond;
paramSave.ioc_cf_q_constraint = param.ioc_cf_q_constraint;
paramSave.ioc_cf_dq_constraint = param.ioc_cf_dq_constraint;
paramSave.ioc_knot_locations = param.ioc_knot_locations;

paramSave.dt_full = param.dt_full;
% paramSave.jump.takeoffFrame = param.jump.takeoffFrame;
% paramSave.jump.landFrame = param.jump.landFrame;
% paramSave.jump.locationLand = param.jump.locationLand;
% paramSave.jump.grade = param.jump.grade;
% paramSave.jump.modelLinks = param.jump.modelLinks;
% paramSave.jump.world2base = param.jump.world2base;
% paramSave.jump.bad2d = param.jump.bad2d;

saveVar.param = paramSave;