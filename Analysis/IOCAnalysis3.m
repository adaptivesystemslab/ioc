function IOCAnalysis()
    basePath = 'D:\aslab_gitlab\expressive-ioc\Data\IOC\plot_20200319162004\';

    currBasePathDir = dir(basePath);
    for j = 1:length(currBasePathDir)
        currSubPath = fullfile(basePath, currBasePathDir(j).name);
        if strcmpi(currSubPath(end), '.')
            continue; % it's . or ..
        end

        % load and plot stuff
        loadAndPlotStuff(currSubPath);
    end
end

function loadAndPlotStuff(filepath)
    load(filepath);
    
end