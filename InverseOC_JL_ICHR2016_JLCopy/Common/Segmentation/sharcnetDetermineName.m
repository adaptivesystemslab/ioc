function export = sharcnetDetermineName(targetFileName)

%     scriptPath = dbstack;
%     targetFileName = scriptPath(end).name;
    currScriptName = strsplit('z', targetFileName);
    exportPathSuffix = currScriptName{2};

    if length(currScriptName) > 2
        pxSet = str2num(currScriptName{3});
    else
        pxSet = 0;
    end
    
    export.exportPathSuffix = exportPathSuffix;
    export.pxSet = pxSet;
end