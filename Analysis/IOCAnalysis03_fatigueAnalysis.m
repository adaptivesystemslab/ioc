function IOCAnalysis()
    setPaths();
%     nowstr = datestr(now, 'yyyymmddHHMMSS');
    nowstr = 'test';
      
    basePath = 'D:\results\fatigue_ioc02_weightsAssembled\plot_20200324120157\';
    searchString = 'mat_dataInd_*.mat';
    filepathSegments = 'ManualSeg.xlsx';
    outputPath = ['D:\results\fatigue_ioc03_weightsPattern\plot_' nowstr '\'];
    checkMkdir(outputPath);
    
    currBasePathDir = dir([basePath searchString]);
    for j = 1:length(currBasePathDir)
        currFileName = currBasePathDir(j).name;
        
        if strcmpi(currFileName(end), '.')
            continue; % it's . or ..
        end

        % load and plot stuff
        filepathCurrDataInd = fullfile(basePath, currFileName);
        filepathCurrWeiCum = strrep(filepathCurrDataInd, 'mat_dataInd', 'mat_weiCum');
        filepathCurrWeiInd = strrep(filepathCurrDataInd, 'mat_dataInd', 'mat_weiInd');
        calculateMetrics(filepathCurrDataInd, filepathCurrWeiCum, filepathCurrWeiInd, filepathSegments, outputPath);
    end
end

function loadAndPlotStuff(filepath)
    load(filepath);
    
end

function calculateMetrics(filepathCurrDataInd, filepathCurrWeiCum, filepathCurrWeiInd, filepathSegments, outputPath)
    load(filepathCurrDataInd);
    load(filepathCurrWeiCum);
    
    % load seg info
    trialInfo = matData.trialInfo;
    [segData, segOnlyDataTable, restOnlyDataTable] = loadSegmentInfo(filepathSegments, trialInfo);
    segmentInfo = segData;
    
    % replace trialInfo's load path
    oldLoadPath = trialInfo.path;
    newLoadPath = strrep(oldLoadPath, 'D:\aslab_svn\data_IK\FullBody_IIT_2017', 'D:\results\fatigue_ioc00_ik');
    trialInfo.path = newLoadPath;
    
    model = getModel(trialInfo); % now need to run FK/FD to get dq and tau
    traj = model.loadData(trialInfo);
    allJointNames = {model.model.joints.name};
    
    for i = 1:size(traj.q, 2)
        featureLabel = ['q_' num2str(i) '_' allJointNames{i}];
        featureCalc1(featureLabel, traj.trajT, traj.q(:, i), segmentInfo, outputPath, ['feat_q_' num2str(i)]);
    end
    for i = 1:size(traj.dq, 2)
        featureLabel = ['dq_' num2str(i) '_' allJointNames{i}];
        featureCalc1(featureLabel, traj.trajT, traj.dq(:, i), segmentInfo, outputPath, ['feat_dq_' num2str(i)]);
    end
    for i = 1:size(traj.tau, 2)
        featureLabel = ['tau_' num2str(i) '_' allJointNames{i}];
        featureCalc1(featureLabel, traj.trajT, traj.tau(:, i), segmentInfo, outputPath, ['feat_tau_' num2str(i)]);
    end
    
    for i = 1:size(matSave.weights, 2)
        featureLabel = ['weights_' matData.featureLabels{i}];
        featureCalc1(featureLabel, matSave.t, matSave.weights(:, i), segmentInfo, outputPath, ['feat_weights_' featureLabel]);
    end
end

function stats = featureCalc1(name, t, feature, segData, outputPath, figFile)
%     segmentInfo = struct2table(segData);
%     mask = ismember(segDataTable{:,'state'}, 'Seg');
%     segmentInfo = segData(mask);

    ind_rest = 0;
    ind_seg = 0;

    % calculate metrics
     for i = 1:length(segData)
        currStartTime = segData(i).timeStart;
        currEndTime = segData(i).timeEnd;
        currState = segData(i).state{1};
        [~, currStartInd] = findClosestValue(currStartTime, t);
        [~, currEndInd] = findClosestValue(currEndTime, t);
        
        if currStartInd < 0
            currStartInd = 1;
        end
        
        if currEndInd > length(feature)
            currEndInd = length(feature);
        end
        
        currFeature = feature(currStartInd:currEndInd, :);
         
        segmentStats.t = (currStartTime+currEndTime)/2;
%         segmentStats.t = i;
        
        segmentStats.mean = mean(currFeature);
        segmentStats.stddev = std(currFeature);
        segmentStats.rms = rms(currFeature);
        segmentStats.skewness = skewness(currFeature);
        segmentStats.kurtosis = kurtosis(currFeature);
        segmentStats.hurst = genhurst(currFeature);
        segmentStats.entropy = wentropy(currFeature,'shannon');
        
        [activity, mobility, complexity] = hjorthParam(currFeature);
        segmentStats.hjorth_activity = activity;
        segmentStats.hjorth_mobility = mobility;
        segmentStats.hjorth_complexity = complexity;
        
        segmentStats.peakFreq = peakFreq(currFeature);
        
        switch currState
            case 'Seg'
                ind_seg = ind_seg + 1;
                segmentStatsSeg(ind_seg) = segmentStats;
                
            case 'Rest'
                ind_rest = ind_rest + 1;
                segmentStatsRest(ind_rest) = segmentStats;
        end
     end
    
     stats.name = name;
     stats.t = t;
     stats.feature = feature;
     stats.segData = segData;
     stats.segmentStatsSeg = segmentStatsSeg;
     stats.segmentStatsRest = segmentStatsRest;
     
     h = plotData(stats);
     figPath = [outputPath, figFile];
     saveas(h, figPath, 'png');
     saveas(h, figPath, 'fig');
     close(h);
end

function h = plotData(stats)
    h = figure('Position', [-1919 69 1920 964.8000]);
    
    segDataTable = struct2table(stats.segData);
    mask = ismember(segDataTable{:,'state'}, 'Seg');
    segOnlyDataTable = stats.segData(mask);
   
    mask = ismember(segDataTable{:,'state'}, 'Rest');
    restOnlyDataTable = stats.segData(mask);

    title(stats.name);
    
    featureParams = {'mean', 'stddev', 'rms', 'skewness', ...
        'kurtosis', 'hurst', 'entropy', 'hjorth_activity', 'hjorth_mobility', 'hjorth_complexity', 'peakFreq'};
    
    for i = 1:length(featureParams)
        ax(i+1) = subplot(4, 3, i+1);
        featureParam = featureParams{i};
        
        segT = [stats.segmentStatsSeg.t];
        segRaw = [stats.segmentStatsSeg.(featureParam)];
        segNorm = segRaw / max(abs(segRaw));
        [b_seg, Rsq2_seg] = linearFit(segT', segNorm');
        
        restT = [stats.segmentStatsRest.t];
        restRaw = [stats.segmentStatsRest.(featureParam)];
        restNorm = restRaw / max(abs(restRaw));
        [b_rest, Rsq2_rest] = linearFit(restT', restNorm');
       
%         yCalc2 = X*b;
%         plot(x,yCalc2,'--')
%         legend('Data','Slope','Slope & Intercept','Location','best');
        
        yyaxis left
        plAx(1) = plot(segT, segNorm, '-o', 'MarkerSize', 16); hold on
        ylabel(['R2_Seg = ' num2str(Rsq2_seg, '%0.2f')]);
        
        yyaxis right
        plAx(2) = plot(restT,	restNorm, '-o', 'MarkerSize', 16);
        ylabel(['R2_Res = ' num2str(Rsq2_rest, '%0.2f')]);
        
        R_thres = 0.7;
        if Rsq2_seg > R_thres && Rsq2_rest > R_thres
            title([featureParam, ' (R2_Res = ' num2str(Rsq2_rest, '%0.2f') ', R2_Seg = ' num2str(Rsq2_seg, '%0.2f') ')']);
        elseif Rsq2_seg > R_thres
            title([featureParam, ' (R2_Seg = ' num2str(Rsq2_seg, '%0.2f') ')']);
        elseif Rsq2_rest > R_thres
            title([featureParam, ' (R2_Res = ' num2str(Rsq2_rest, '%0.2f') ')']);
        else
            title(featureParam);
        end
        
%         ylim([-1 1]);
    end
    
    ax(1) = subplot(4, 3, 1);
    plot(stats.t, stats.feature);
    plotBoxes(segOnlyDataTable, plAx(1).Color);
    plotBoxes(restOnlyDataTable, plAx(2).Color);
    
    linkaxes(ax, 'x');
end

function [b, Rsq2] = linearFit(x, y)
    X = [ones(length(x),1) x];
    b = X\y;
    yCalc2 = X*b;
    Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
end

function frqs = peakFreq(signal)
%     FTsignal = fft(signal - mean(signal))/length(signal);
%     [maxpeak, maxpeakindes] = max(abs(FTsignal)*2);

    F=signal;                % Data Channel
    Ts = 0.01;                  % Sampling Interval (s)
    Fs = 1/Ts;                  % Sampling Frequency (Hz)
    Fn = Fs/2;                  % Nyquist Frequency
    F(isnan(F))=[];             % Eliminate �NaN� Values First
    LF = size(F,1);             % Length of Data Vector
    T = linspace(0,1,LF)*Ts;    % Create Time Vector
 
%     figure(1)                   % Plot Data
%     plot(T, F)
%     grid
    FF = fft(F)/LF;             % Fourier Series of Data, Freq Vector
    Fv = linspace(0,1,fix(LF/2)+1)*Fn;
    Iv = 1:length(Fv);          % Index Vector
   
%     figure(2)                   % Plot FFT
%     plot(Fv, abs(FF(Iv)))
%     grid
%     xlabel('Frequency (Hz)')
%     ylabel('Amplitude')
%     axis([0  1500    ylim])
    x = abs(FF);
  
%     figure()
    [pks, locs]=findpeaks(x(Iv));
%     plot(Fv(locs), pks, 'or')
    
if ~isempty(locs)
    frqs = Fv(locs(1));
else
    frqs = 0;
end
%     hold on;
%     plot(x)
%     [pks, locs]=findpeaks(x, 'MinPeakDistance',50, 'minpeakheight',0.002);
end

function [activity, mobility, complexity] = hjorthParam(currFeature)
        %% Hjorth parameters
    % Activity
    activity = var(currFeature);
    
    % Mobility
    mobility = std(diff(currFeature))./std(currFeature);
    
    % Complexity
    complexity = std(diff(diff(currFeature)))./std(diff(currFeature))./mobility;
end

function feat = calculate_features(f)
% CALCULATE_FEATURES computes a set of features that will be
%    later on used to train a LASSO GLM model.
%
%    The set of features include: spectral entropy and Shannons
%    entropy (MacKay, 2003) at six frequency bands: delta (0.14 Hz),
%    theta (48 Hz), alpha (812 Hz), beta (1230 Hz), low-gamma (3070 Hz)
%    and high gamma (70180 Hz), and Shannons entropy in dyadic (between
%    0.00167 and 109 Hz spaced by factors of 2n) frequency bands, the
%    spectral edge at 50% power below 40 Hz, spectral correlation between
%    channels in dyadic frequency bands, the time series correlation matrix
%    and its eigenvalues, fractal dimensions, Hjorth activity, mobility and
%    complexity parameters (Hjorth, 1970), and the statistical skewness and
%    kurtosis of the distribution of time series values.
%
%    Input:
%       f: data structure that includes data and metadata information
%          e.g. f = load('1_1_1.mat');
%   
%    Output:
%       feat: a vector of the abovementioned features
%
%    Reference:
%       http://brain.oxfordjournals.org/content/early/2016/03/29/brain.aww045
%       https://github.com/drewabbot/kaggle-seizure-prediction
%
%  Copyright 2016 The MathWorks, Inc.
%

fName = fieldnames(f);
fs = f.(fName{1}).iEEGsamplingRate;     % Sampling Freq
eegData = f.(fName{1}).data;            % EEG data matrix
[nt,nc] = size(eegData);

subsampLen = floor(fs * 60);            % Num samples in 1 min window
numSamps = floor(nt / subsampLen);      % Num of 1-min samples
sampIdx = 1:subsampLen:numSamps*subsampLen;

feat = []; % Feature vector

for l=2:numSamps
    %% Sample 1-min window
    epoch = eegData(sampIdx(l-1):sampIdx(l),:);
    
    %% Compute Shannon's entropy, spectral edge and correlation matrix
    % Find the power spectrum at each frequency bands
    D = abs(fft(eegData));                  % take FFT of each channel
    D(1,:) = 0;                             % set DC component to 0
    D = bsxfun(@rdivide,D,sum(D));          % normalize each channel
    lvl = [0.1 4 8 12 30 70 180];           % frequency levels in Hz
    lseg = round(nt/fs*lvl)+1;              % segments corresponding to frequency bands
    
    dspect = zeros(length(lvl)-1,nc);
    for n=1:length(lvl)-1
        dspect(n,:) = 2*sum(D(lseg(n):lseg(n+1),:));
    end
    
    % Find the Shannon's entropy
    spentropy = -sum(dspect.*log(dspect));
    
    % Find the spectral edge frequency
    sfreq = fs;
    tfreq = 40;
    ppow = 0.5;
    
    topfreq = round(nt/sfreq*tfreq)+1;
    A = cumsum(D(1:topfreq,:));
    B = bsxfun(@minus,A,max(A)*ppow);
    [~,spedge] = min(abs(B));
    spedge = (spedge-1)/(topfreq-1)*tfreq;
    
    % Calculate correlation matrix and its eigenvalues (b/w channels)
    type_corr = 'Pearson';
    C = corr(dspect,'type',type_corr);
    C(isnan(C)) = 0;    % make NaN become 0
    C(isinf(C)) = 0;    % make Inf become 0
    lxchannels = real(sort(eig(C)));
    
    % Calculate correlation matrix and its eigenvalues (b/w freq)
    C = corr(dspect.','type',type_corr);
    C(isnan(C)) = 0;    % make NaN become 0
    C(isinf(C)) = 0;    % make Inf become 0
    lxfreqbands = real(sort(eig(C)));
    
    %% Spectral entropy for dyadic bands
    % Find number of dyadic levels
    ldat = floor(nt/2);
    no_levels = floor(log2(ldat));
    seg = floor(ldat/2^(no_levels-1));
    
    % Find the power spectrum at each dyadic level
    dspect = zeros(no_levels,nc);
    for n=no_levels:-1:1
        dspect(n,:) = 2*sum(D(floor(ldat/2)+1:ldat,:));
        ldat = floor(ldat/2);
    end
    
    % Find the Shannon's entropy
    spentropyDyd = -sum(dspect.*log(dspect));
    
    % Find correlation between channels
    C = corr(dspect,'type',type_corr);
    C(isnan(C)) = 0;    % make NaN become 0
    C(isinf(C)) = 0;    % make Inf become 0
    lxchannelsDyd = sort(eig(C));
    
    %% Fractal dimensions
    no_channels = nc;
    fd = zeros(3,no_channels);
    
    for n=1:no_channels
        fd(:,n) = wfbmesti(eegData(:,n));
    end
    
    %% Hjorth parameters
    % Activity
    activity = var(eegData);
    
    % Mobility
    mobility = std(diff(eegData))./std(eegData);
    
    % Complexity
    complexity = std(diff(diff(eegData)))./std(diff(eegData))./mobility;
    
    %% Statistical properties
    % Skewness
    skew = skewness(eegData);
    
    % Kurtosis
    kurt = kurtosis(eegData);
    
    %% Compile all the features
    feat = [feat spentropy(:)' spedge(:)' lxchannels(:)' lxfreqbands(:)' spentropyDyd(:)' ...
        lxchannelsDyd(:)' fd(:)' activity(:)' mobility(:)' complexity(:)' skew(:)' kurt(:)'];
end
end

% function perSegmentMetrics()
%     for i = 1:length(segmentInfo)
%         startTime = segmentInfo(i).
%     end
% end