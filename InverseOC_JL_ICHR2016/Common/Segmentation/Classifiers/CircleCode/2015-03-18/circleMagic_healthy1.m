function circleMagic_healthy1

% set up specStruct, which determines the specifics of the data you want to
% load
specStruct.datasetName = 'healthy1'; % the dataset name: ie 'healthy1', 'healthy2'
specStruct.patient = [3]; % type in px numbers here: ie [1 13 14 17 18]. Leave blank if want all px
specStruct.session = [1]; % type in the session number here: ie [1 4 6]. Leave blank if want all sessions

specStruct.exerciseAcceptPrefix = {'KEFO_SIT_SLO'}; % exercise prefix that you want. This will load only exercises that start with 'KEFO_SIT'
specStruct.exerciseAcceptSuffix = {}; % exercise suffix that you want. This will load only exercises that end with 'SLO'
specStruct.exerciseRejectPrefix = {};  % exercise prefix that you don't want
specStruct.exerciseRejectSuffix = {};  % exercise suffix that you don't wantwant.

% options.headerMocap = []; % so load mocap and imu data
% options.headerImu = [];
% options.headerEkf = []; % but don't load ekf and seg data
options.headerSegManual = [];
options.headerSegCrop = [];
options.headerSegAlg = [];

pathToRawData = callDatasetBasePath(specStruct); % update the contents of this function with local specifications
fileStack = loadPatientFilepaths(pathToRawData, specStruct);

    % load data
%     exerciseStruct = cell(size(fileStack));
    for ind_fileStack = 1
        exerciseStruct{ind_fileStack} = exerciseDataHandle(fileStack{ind_fileStack}.filePath, options);
        exerciseStruct{ind_fileStack}.load;
    end

% now start taking things apart
X = exerciseStruct{1}.dataEkf.time - exerciseStruct{1}.dataEkf.time(1);
Y = exerciseStruct{1}.dataEkf.Q(:, 4);

[myFit, Y_sin, Y_cos] = nonLinearModelFit(X, Y);
% myFit = fourierSeriesFit(X, Y);

%% Generate a plot
figure
scatter(X, Y, 'r')
hold on
plot(X(300:500), myFit.Fitted)
hold off

figure
plot(Y_sin, Y_cos);
% plot(Y_cos, Y_sin);

function [myFit, Y_sin, Y_cos] = nonLinearModelFit(X, Y)
    % want to fit to a sum of y = b0 + Ai*sin(w *t + phii) +  Ai*sin(w *t + phii)
    %                             b0 + b1*sin(b2*t + b3)
    
    % scale down the X and Y
    indRange = 300:500;
    X = X(indRange);
    Y = Y(indRange);
    
    w = [0.1:0.1:4];      % b2x
    A = ones(size(w));    % b1x
    phi = zeros(size(w)); % b3x

    fittingStr = 'y ~ b0 + ';
    fittingStr2 = 'y ~ b0 + ';
    fittingInput = [-1.5];
    for ind = 1:length(A)
        n2si = num2str(ind);
        fittingStr = [fittingStr 'b1' n2si '*sin(' num2str(w(ind)) '*x1+b3' n2si ') + ']; % fixed freq intervals
        fittingStr2 = [fittingStr2 'b1' n2si '*cos(' num2str(w(ind)) '*x1+b3' n2si ') + '];
        fittingInput = [fittingInput A(ind) phi(ind)];
%         fittingStr = [fittingStr 'b1' n2si '*sin(b2' n2si '*x1+b3' n2si ') + '];
%         fittingInput = [fittingInput A(ind) w(ind) phi(ind)];
    end
    
    fittingStr = fittingStr(1:end-2);
    fittingStr2 = fittingStr2(1:end-2);
    
    myFit = NonLinearModel.fit(X,Y, fittingStr, fittingInput)
    
%     now sub the equations back and generate the y
    Y_sin = ones(size(X))*myFit.Coefficients{'b0', 'Estimate'};
    Y_cos = ones(size(X))*myFit.Coefficients{'b0', 'Estimate'};
    for ind = 1:length(A)
        n2si = num2str(ind);
        b1 = myFit.Coefficients{['b1' n2si], 'Estimate'};
        b3 = myFit.Coefficients{['b3' n2si], 'Estimate'};
        Y_sin = Y_sin + b1*sin(w(ind)*X+b3);
        Y_cos = Y_cos + b1*cos(w(ind)*X+b3);
    end

    
function myFit = fourierSeriesFit(X, Y)
    f2 = fit(X,Y,'fourier8')
    plot(f2,X,Y)
