feature_load_names{1} = 't';
feature_load_names{2} = 'q';
feature_load_names{3} = 'dq';
feature_load_names{4} = 'ddq';
feature_load_names{5} = 'dddq';
feature_load_names{6} = 't_cf_q_const';
feature_load_names{7} = 't_cf_dq_const';
feature_load_names{8} = 't_cf_ddq_const';
feature_load_names{9} = 't_knot';
feature_load_names{10} = 't_sp_dq_const';
feature_load_names{11} = 'y_cf_q_const';
feature_load_names{12} = 'y_cf_dq_const';
feature_load_names{13} = 'y_cf_ddq_const';
feature_load_names{14} = 'y_knot';
feature_load_names{15} = 'y_sp_dq_const';

feature_full = feature_collate(traj_load_lengthAdj, feature_load_names);

%     y_opt_mat = vec2mat_multi(y_opt_vec, param.matVec_struct); % convert the vector back into matrix form for reading
%     splineFit.q = calc_direct_spline(y_opt_mat, param); % calc the spline from knot points (ie the variables to opt for)
%     
%     feature_use = calc_features(splineFit, feature_win, param);

plotKnots = 1;
plotSegments = 0;

clc
exitStr = '';
for i = 1:length(output_optsim)
    output_optsim(i)
    exitStr = [exitStr ', ' num2str(output_optsim(i).exitflag)];
end

 figure;
        ax(1) = subplot(221);
        plot(feature_full.t, feature_full.q'); title('q');
        if plotSegments
            hold on
            for i = [0 2 4 6 8]
                plot([i i], ylim, 'k');
            end
        end        
        title(['exitflag: ' exitStr]);
        
        ax(2) = subplot(222);
        plot(feature_full.t, feature_full.dq'); title('dq');
%         hold on; plot(dq', 'r');
        ax(3) = subplot(223);
        plot(feature_full.t, feature_full.ddq'); title('ddq');

        if plotKnots
            hold on
            for i = 1:length(feature_full.t_knot)
                plot([feature_full.t_knot(i) feature_full.t_knot(i)], ylim, 'k');
            end
        end
        
        ax(4) = subplot(224);
        plot(feature_full.t, feature_full.dddq'); title('dddq');
        
                if plotKnots
            hold on
            for i = 1:length(feature_full.t_knot)
                plot([feature_full.t_knot(i) feature_full.t_knot(i)], ylim, 'k');
            end
                end
        
                linkaxes(ax, 'x');