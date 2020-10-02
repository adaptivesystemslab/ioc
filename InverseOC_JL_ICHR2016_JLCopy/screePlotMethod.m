function [dof, eig_vals, eig_perc] = screePlotMethod(array, threshold)
    
    eig_vals = svd(array);
    eig_perc = eig_vals/sum(eig_vals);
    
    sum_eig = 0;
    index = 0;
    while sum_eig < threshold && index < length(eig_vals)
        index = index + 1;
        sum_eig = sum_eig + eig_perc(index);
    end
    dof = index;
    
    eig_perc = eig_perc*100;
