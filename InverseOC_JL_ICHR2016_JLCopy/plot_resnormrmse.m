function h = plot_resnormrmse(output_inverse, minRmseIndArray, indToUse_window, temp)
% plot rmse to resnorm

    resnorm = [];
    rmse = [];
    
%     output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)))
    windowCount = size(indToUse_window, 1);
    
    for ind_windowCount = 1:windowCount
        output_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount));
        
        resnorm = [resnorm output_curr.resnorm];
        rmse = [rmse output_curr.rmse];
%         rmse = [rmse temp(ind_windowCount)];
    end
    
    h = figure;
    plot(rmse, resnorm, 'rx');
    xlabel('rmse');
    ylabel('resnorm');
    hold on
    
    % filter it a bit
    indToRemove1 = find(rmse > 10);
    indToRemove2 = find(resnorm > 10);
    indToRemove = unique([indToRemove1 indToRemove2]);
    indTotal = 1:length(rmse);
    indToKeep = setxor(indTotal, indToRemove);
    
    rmse = rmse(indToKeep) ;
    resnorm = resnorm(indToKeep);
    
    p = polyfit(rmse, resnorm, 1);
    f = polyval(p,rmse);
    
    plot(rmse, resnorm, 'bo');
    plot(rmse, f, 'b-');
    
    x = deg2rad(5);
    f = polyval(p,x);
    plot([x x], ylim, 'k-');
    title(['Threshold for RMSE = ' num2str(x) ' is RESNORM ' num2str(f)]);
end