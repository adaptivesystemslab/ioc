function fileConverter(directory, sourceFileType, targetFileType)
    % converts files of one type into another
    
    if ~exist('directory', 'var')
        directory = 'D:\MATLABResults\EKF5DOF\2012-11-30_delwhendone\EKF5\';
        targetDir = 'D:\MATLABResults\EKF5DOF\2012-11-30_delwhendone\EKF5\2238\';
        sourceFileType = 'fig';
        targetFileType = 'eps';
        targetFileTypeSpec = 'epsc';
    end
    
    if ~exist(targetDir)
        mkdir(targetDir); 
    end
    
    baseFolderDir = dir(directory);
    for i = 1:length(baseFolderDir)
        currFolderStruct = baseFolderDir(i);
        fullFilePath = [directory currFolderStruct.name];
        
        % check the currently 'active' folder
        if strcmp(currFolderStruct.name(1), '.')
            % if the dir result starts with a period...not wanted
            continue
        elseif currFolderStruct.isdir == 1
            continue
        elseif ~strcmp(currFolderStruct.name(end-2:end), sourceFileType)
            % not the type of file we're looking for
            continue
        end
        
        h = hgload(fullFilePath);
        
        % if the labels needs to be modified
%         xlim([6 20]);
        title('');
%         xlabel('Time [s]');
%         ylabel('End-effector position [m]');
        
N = 24;
set(gca,'fontsize',N)
set(findall(gcf,'type','text'),'fontSize',N)

% extract information for this instance (ekf things)
if ~isempty(str2num(currFolderStruct.name(2)))
    % this is 10
    mvtNumber = currFolderStruct.name(1:2);
    setNumber = currFolderStruct.name(8);
else
        mvtNumber = currFolderStruct.name(1);
    setNumber = currFolderStruct.name(7);
end
   position = currFolderStruct.name(end-9:end-5);
    figType = currFolderStruct.name(end-4);
    
    xlim([6 20])
    
    switch figType
        case '1'
            figTypeStr = 'accel';
        case '2'
            figTypeStr = 'gyro';
            ylim([-2.2 2.2])
            
        case '3'
            figTypeStr = 'joint';
        case '4'
            figTypeStr = 'ef';
    end
    
    
    
output = ['5_' mvtNumber '_slow' setNumber '_ekf_' position '_' figTypeStr];

        saveas(h, [targetDir output '.' targetFileType], targetFileTypeSpec);
        
% hmm things
% xlim([22 38])
% 
%         saveas(h, [targetDir currFolderStruct.name(1:end-4) '.' targetFileType], targetFileTypeSpec);



        close(h);
    end
end