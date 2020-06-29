% loading example for exerciseLoader
if exist('specStruct', 'var')
    clear specStruct options exerciseStruct
end

% set up specStruct, which determines the specifics of the data you want to
% load
specStruct.datasetName = 'healthy1'; % the dataset name: ie 'healthy1', 'healthy2'
specStruct.patient = [3]; % type in px numbers here: ie [1 13 14 17 18]. Leave blank if want all px
specStruct.session = [1]; % type in the session number here: ie [1 4 6]. Leave blank if want all sessions

specStruct.exerciseAcceptPrefix = {'KEFO_SIT_SLO1'}; % exercise prefix that you want. This will load only exercises that start with 'KEFO_SIT'
% specStruct.exerciseAcceptSuffix = {}; % exercise suffix that you want. This will load only exercises that end with 'SLO'
% specStruct.exerciseRejectPrefix = {};  % exercise prefix that you don't want
% specStruct.exerciseRejectSuffix = {};  % exercise suffix that you don't wantwant.

options.headerMocap = []; % so load mocap and imu data
options.headerImu = [];
% options.headerEkf = []; % but don't load ekf and seg data
options.headerSegManual = [];
options.headerSegCrop = [];
options.headerSegAlg = [];

exerciseStruct = exerciseLoader(specStruct, options);

exerciseStruct{1}
exerciseStruct{1}.plot % print to console for demo
exerciseStruct{1}.cropData
% exerciseStruct{1}.plot % print to console for demo