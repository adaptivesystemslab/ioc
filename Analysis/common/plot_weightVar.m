function plot_weightVar(weights_all, weights_all_var, t,q,faceColours, outputPathFig_all_variance, weightLabels)
 h = figure('Position', [488 342 560*2 420*2]);
    for i = 1:size(weights_all, 2)+1
        switch size(weights_all, 2)
            case 3
                ax(i) = subplot(2,2,i);
            case 4
                ax(i) = subplot(2,3,i);
            case 8
                ax(i) = subplot(3,3,i);
        end
        
        switch i
            case 1
                plot(t, q);
                title('q');
            otherwise
                ind_weight = i-1;
%                 boundedline(t,weights_all(:,ind_weight), weights_all_var(:,ind_weight), 'cmap', faceColours(ind_weight, :));
                boundedline(t,weights_all(:,ind_weight), weights_all_var(:,ind_weight), '-b');
                title(weightLabels{ind_weight});
                ylim([-0.25 1.25]);
        end
    end
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig_all_variance, 'fig');
    saveas(h, outputPathFig_all_variance, 'png');
    close(h);
end