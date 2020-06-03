% loading example for exerciseLoader
if exist('specStruct', 'var')
    clear specStruct options
end

% set up specStruct, which determines the specifics of the data you want to
% load
specStruct.datasetName = 'healthy1'; % the dataset name: ie 'healthy1', 'healthy2'
specStruct.patient = [1]; % type in px numbers here: ie [1 13 14 17 18]. Leave blank if want all px
specStruct.session = []; % type in the session number here: ie [1 4 6]. Leave blank if want all sessions

% the following accepts/rejects specific exercise types that it comes
% across. The length of the string it checks against will be the same as
% the length as the first entry in each cell array, so please be consist in
% string length. Leave these arrays blank if not needed. Different 
% combinations of this section will achieve the same results. For example, 
% the following two lines can also be replaced with:
%     exerciseAcceptPrefix = {'KEFO_SIT_SLO'};
% and it would achieve the same thing
specStruct.exerciseAcceptPrefix = {'KEFO_SIT_SLO1'}; % exercise prefix that you want. This will load only exercises that start with 'KEFO_SIT'
specStruct.exerciseAcceptSuffix = {}; % exercise suffix that you want. This will load only exercises that end with 'SLO'
specStruct.exerciseRejectPrefix = {};  % exercise prefix that you don't want
specStruct.exerciseRejectSuffix = {};  % exercise suffix that you don't wantwant.

% by default, the exerciseLoader looks at default locations and attempts to
% load all available data that is associated to a given
% patient/session/exercise combination. If you don't want to load the
% default paths to EKF/IMU/mocap etc, you can pass in specific paths for
% the exerciseDataHandle to take (which, however, means you will have to
% write your own exerciseLoader function). 

% However, exerciseLoader can accept arguments to stop loading a given data
% type of data. If you don't want to load a given set of data, specify the 
% type and pass in an empty array. For example, the following will ensure 
% that efk, and all segmentation data will not be loaded.

% Adjust the options array appropriately to suit your needs. If no options
% file is passed in, all available data is used loaded

options.headerMocap = []; % so load mocap and imu data
options.headerImu = [];
% options.headerEkf = []; % but don't load ekf and seg data
% options.headerSegManual = [];
options.headerSegCrop = [];
options.headerSegAlg = [];

% example for exerciseLoader
exerciseStruct = exerciseLoader(specStruct, options);
exerciseStruct{1}
exerciseStruct{1}.plot % print to console for demo
exerciseStruct{1}.cropData % apply the cropping from seg crop, if seg crop exists