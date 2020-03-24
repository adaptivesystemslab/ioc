function IOCAnalysis()
    setPaths();
    
    basePath = 'D:\results\fatigue_ioc02_weightsAssembled\plot_20200324120157\';
    searchString = 'mat_dataInd_*.mat';
    filepathSegments = 'ManualSeg.xlsx';
    outputPath = 'D:\results\fatigue_ioc02_weightsAssembled';
    
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
        calculateMetrics(filepathCurrDataInd, filepathCurrWeiCum, filepathCurrWeiInd, filepathSegments);
    end
end

function loadAndPlotStuff(filepath)
    load(filepath);
    
end

function calculateMetrics(filepathCurrDataInd, filepathCurrWeiCum, filepathCurrWeiInd, filepathSegments)
    load(filepathCurrDataInd);
    load(filepathCurrWeiCum);
    
    % load seg info
    trialInfo = matData.trialInfo;
    [segData, segOnlyDataTable, restOnlyDataTable] = loadSegmentInfo(filepathSegments, trialInfo);
    
    segmentInfo = segOnlyDataTable;
    t = matData.t;
    feature = matData.q(:, 1);
    featureCalc1('q', t, feature, segmentInfo);
end

function stats = featureCalc1(name, t, feature, segmentInfo)
    % calculate metrics
     for i = 1:length(segmentInfo)
        currStartTime = segmentInfo(i).timeStart;
        currEndTime = segmentInfo(i).timeEnd;
        [~, currStartInd] = findClosestValue(currStartTime, t);
        [~, currEndInd] = findClosestValue(currEndTime, t);
        currFeature = feature(currStartInd:currEndInd, :);
         
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
        
%         segmentStats.freq = peakFreq(currFeature);
        segmentStatsOut(i) = segmentStats;
     end
    
     stats.name = name;
     stats.t = t;
     stats.feature = feature;
     stats.segment = segmentInfo;
end

function plotData(stats)

end

% function frqs = peakFreq(currFeature)
% Fs = 100;
% L = length(currFeature);             % Length of signal
% 
% X = currFeature;
% Y = fft(X);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
%     freq = fft(abs(currFeature));
%     [pks,frqs] = findpeaks(abs(currFeature),freq);
% end

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