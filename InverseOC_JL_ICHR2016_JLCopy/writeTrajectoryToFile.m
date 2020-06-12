function writeTrajectoryToFile(exportPath, currInstName, cost_function_names, time, avgWeight, avgRatio)
    filename = fullfile(exportPath, [currInstName '_Trajectory.txt']);
    if ~exist(filename, 'file') 
        % if doesn't exist, write header  
        header = ["time"];
        outputStr = ['%s\t'];
        
        for i = 1:length(cost_function_names)
            header = [header,string(['Weight_' cost_function_names{i}])];
            outputStr = [outputStr '%s\t'];
        end
        
        for i = 1:length(cost_function_names)
            header = [header, string(['Ratio_' cost_function_names{i}])];
            outputStr = [outputStr '%s\t'];
        end
        outputStr = [outputStr '\n'];
    else
        header = '';
    end
    
    fid = fopen(filename, 'wt');
    fprintf(fid, outputStr, header);
    fclose(fid);
    
    keepInds = find(time > 0);
    t_obs = time(keepInds);
    weight_obs = avgWeight(keepInds, :);
    ratio_obs = avgRatio(keepInds, :);  
    
    data = [t_obs', weight_obs, ratio_obs];
    dlmwrite(filename,data,'delimiter','\t','precision',['%10.',num2str(12),'f'],'newline', 'unix', '-append');
end